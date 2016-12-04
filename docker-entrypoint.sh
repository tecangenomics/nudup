#!/bin/bash

echo "Python version"
python --version;

echo "Docker build version"
cat /version;

cd /opt/nudup

echo "Requires test_data to complete test"
python test.py
