#!/usr/bin/env bash
# Packages the release build as dist/lrsqueeze-VERSION-ARCH-OS.tar.gz, with the binary, README,
# LICENSE and CHANGELOG, and its SHA-256 checksum beside it.
# usage: scripts/package-release.sh VERSION   (run from the repository root after cargo build --release)
set -euo pipefail
version=$1
case "$(uname -s)" in
    Linux) os=linux ;;
    Darwin) os=macos ;;
    *) echo "unsupported OS $(uname -s)" >&2; exit 1 ;;
esac
arch=$(uname -m)
[ "$arch" = arm64 ] && arch=aarch64
name=lrsqueeze-$version-$arch-$os
rm -rf "dist/$name"
mkdir -p "dist/$name"
cp target/release/lrsqueeze README.md LICENSE CHANGELOG.md "dist/$name/"
tar -C dist -czf "dist/$name.tar.gz" "$name"
if command -v sha256sum > /dev/null; then
    (cd dist && sha256sum "$name.tar.gz" > "$name.tar.gz.sha256")
else
    (cd dist && shasum -a 256 "$name.tar.gz" > "$name.tar.gz.sha256")
fi
echo "dist/$name.tar.gz"
