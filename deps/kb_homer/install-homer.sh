#!/bin/bash

mkdir -p homer
cd homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
cp ./bin/* /kb/deployment/bin
cd ..
