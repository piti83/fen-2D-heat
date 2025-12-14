#!/bin/bash

set +x
echo "Configuring in Debug mode..."

cmake -DCMAKE_BUILD_TYPE=Debug -B build -S .

echo "Configuring done."
