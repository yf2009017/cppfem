#!/bin/sh
cmake . && make && valgrind -q ./cppfem
