//! An ordered parallel map: one thread reads batches, workers transform them, and the calling
//! thread writes the results in input order.

use std::{
    collections::BTreeMap,
    sync::{
        Arc, Mutex,
        mpsc::{self, Receiver, SyncSender},
    },
    thread,
};

use anyhow::ensure;

/// A transformed batch, numbered in input order, or why it failed.
type Done<T> = anyhow::Result<(u64, T)>;

/// Runs `read` on its own thread, `work` on one thread per worker state in `workers`, and `write`
/// on the calling thread, which gets each batch's result in the order `read` sent the batch.
/// `read` sends batches with the function it is given, which returns `false` once the workers
/// have stopped (`read` should then return; the error that stopped them is reported). Returns the
/// worker states, for whatever they counted.
///
/// # Errors
///
/// Fails with the first error from `read`, `work`, or `write`, or if a thread panicked. Every
/// thread has stopped when this returns.
///
/// # Panics
///
/// If `workers` is empty.
pub(crate) fn ordered<B: Send, T: Send, S: Send>(
    read: impl FnOnce(&mut dyn FnMut(B) -> bool) -> anyhow::Result<()> + Send,
    workers: Vec<S>,
    work: impl Fn(&mut S, B) -> anyhow::Result<T> + Sync,
    mut write: impl FnMut(T) -> anyhow::Result<()>,
) -> anyhow::Result<Vec<S>> {
    let threads = workers.len();
    assert!(threads > 0, "the pipeline needs a worker");
    thread::scope(|scope| {
        let (batches, jobs) = mpsc::sync_channel::<(u64, B)>(threads);
        let (results_tx, results) = mpsc::sync_channel(threads);
        let reader = scope.spawn(move || {
            let mut index = 0;
            read(&mut |batch| {
                let sent = batches.send((index, batch)).is_ok();
                index += 1;
                sent
            })
        });
        let jobs = Arc::new(Mutex::new(jobs));
        let work = &work;
        let handles: Vec<_> = workers
            .into_iter()
            .map(|mut state| {
                let (jobs, results) = (Arc::clone(&jobs), results_tx.clone());
                scope.spawn(move || {
                    let _notice = PanicNotice(&results);
                    loop {
                        let Ok((index, batch)) = jobs.lock().expect("lock poisoned").recv() else {
                            break;
                        };
                        let result = work(&mut state, batch).map(|done| (index, done));
                        let failed = result.is_err();
                        if results.send(result).is_err() || failed {
                            break;
                        }
                    }
                    state
                })
            })
            .collect();
        // Only the workers may hold these, so the channels close when they stop.
        drop((jobs, results_tx));

        // Returning from write_in_order drops `results`, which stops the workers and then the
        // reader; joining every thread here turns a panic into an error.
        let written = write_in_order(results, &mut write);
        let read = reader.join();
        let mut states = Vec::with_capacity(handles.len());
        let mut panicked = false;
        for handle in handles {
            match handle.join() {
                Ok(state) => states.push(state),
                Err(_) => panicked = true,
            }
        }
        ensure!(!panicked, "a transform worker panicked");
        written?;
        read.map_err(|_| anyhow::anyhow!("the reader thread panicked"))??;
        Ok(states)
    })
}

/// Tells the writer that a worker panicked, so it stops instead of waiting for the lost batch.
struct PanicNotice<'a, T>(&'a SyncSender<Done<T>>);

impl<T> Drop for PanicNotice<'_, T> {
    fn drop(&mut self) {
        if thread::panicking() {
            // A failed send means the writer has already stopped.
            let _ = self
                .0
                .send(Err(anyhow::anyhow!("a transform worker panicked")));
        }
    }
}

/// Passes transformed batches to `write` in input order. Returning drops `results`, which stops
/// the workers.
fn write_in_order<T>(
    results: Receiver<Done<T>>,
    write: &mut impl FnMut(T) -> anyhow::Result<()>,
) -> anyhow::Result<()> {
    let mut pending = BTreeMap::new();
    let mut next = 0;
    for result in results {
        let (index, done) = result?;
        pending.insert(index, done);
        while let Some(done) = pending.remove(&next) {
            write(done)?;
            next += 1;
        }
    }
    ensure!(pending.is_empty(), "a transform worker stopped early");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Doubles 0..n in batches of one on 4 workers, each counting its batches.
    fn double(
        n: u64,
        work: impl Fn(&mut u64, u64) -> anyhow::Result<u64> + Sync,
        write: impl FnMut(u64) -> anyhow::Result<()>,
    ) -> anyhow::Result<Vec<u64>> {
        ordered(
            |send| {
                for i in 0..n {
                    if !send(i) {
                        break;
                    }
                }
                Ok(())
            },
            vec![0; 4],
            work,
            write,
        )
    }

    #[test]
    fn results_are_written_in_input_order() -> anyhow::Result<()> {
        let mut out = Vec::new();
        let counts = double(
            1000,
            |count, i| {
                *count += 1;
                Ok(2 * i)
            },
            |x| {
                out.push(x);
                Ok(())
            },
        )?;
        assert_eq!(out, (0..1000).map(|i| 2 * i).collect::<Vec<_>>());
        assert_eq!(counts.iter().sum::<u64>(), 1000);
        Ok(())
    }

    #[test]
    fn a_failing_batch_or_write_stops_every_thread() {
        let error = double(
            1000,
            |_, i| {
                ensure!(i != 500, "bad batch {i}");
                Ok(i)
            },
            |_| Ok(()),
        )
        .expect_err("batch 500 fails");
        assert_eq!(error.to_string(), "bad batch 500");

        let mut written = 0;
        let error = double(
            1000,
            |_, i| Ok(i),
            |_| {
                written += 1;
                ensure!(written < 50, "disk full");
                Ok(())
            },
        )
        .expect_err("the 50th write fails");
        assert_eq!(error.to_string(), "disk full");
    }

    #[test]
    fn a_panicking_worker_is_an_error() {
        let error = double(
            1000,
            |_, i| {
                assert!(i != 500, "a bug in a transform");
                Ok(i)
            },
            |_| Ok(()),
        )
        .expect_err("a worker panics");
        assert_eq!(error.to_string(), "a transform worker panicked");
    }

    #[test]
    fn a_failing_reader_is_an_error() {
        let error = ordered(
            |send: &mut dyn FnMut(u64) -> bool| {
                send(0);
                anyhow::bail!("truncated input")
            },
            vec![(); 2],
            |(), i| Ok(i),
            |_| Ok(()),
        )
        .expect_err("the reader fails");
        assert_eq!(error.to_string(), "truncated input");
    }
}
