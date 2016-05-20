#!/bin/sh

# Standard build:
PREFIX=""
:${PREFIX:=$1}
STATIC="-static"

${PREFIX}gcc -Os -Wall -I.. -c ssl_helper.c -o ssl_helper.o
${PREFIX}gcc $STATIC -Wl,--start-group ssl_helper.o -lm ../src/.libs/libwolfssl.a -Wl,--end-group -o ssl_helper
${PREFIX}strip ssl_helper
