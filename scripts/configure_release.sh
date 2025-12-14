#!/bin/bash

set +x
echo "Configuring in Release mode..."

cmake -DCMAKE_BUILD_TYPE=Release -B build -S .

echo "Configuring done."
