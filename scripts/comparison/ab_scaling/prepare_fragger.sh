#!/bin/sh
set -eu

dest=${1:-/tmp/fragger-ab-benchmark}

if [ -e "$dest" ]; then
    echo "destination already exists: $dest" >&2
    exit 2
fi

git clone https://github.com/rhz/graph-rewriting.git "$dest"
git -C "$dest" checkout --detach 0d9595c5d6e8918fc605f7d65778568a74f5e2d4
mkdir -p "$dest/compiler"

for part in compiler library reflect; do
    curl -fL \
      "https://repo.maven.apache.org/maven2/org/scala-lang/scala-$part/2.13.2/scala-$part-2.13.2.jar" \
      -o "$dest/compiler/scala-$part-2.13.2.jar"
done

printf '%s  %s\n' \
  b04da83111decff5b36e87c74d73f5e2231b63ce605f2602e862bae917d7761a \
  "$dest/compiler/scala-compiler-2.13.2.jar" \
  6a1c9846feb8f5d7d33b91e61e056d272215a79747a2b7b5e5ece690bf151d5d \
  "$dest/compiler/scala-library-2.13.2.jar" \
  21a7288874f186d6049e8b14a986011ffd96cb587f206eccf946624f2d476815 \
  "$dest/compiler/scala-reflect-2.13.2.jar" \
  | sha256sum -c -

echo "$dest"
