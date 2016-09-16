#!/bin/sh
yamtbx.shika_backend nproc=0 bl=32xu > ~/shika_backend.log 2>&1 &
cheetah.streaming_client > ~/shika_cheetah.log 2>&1 &
wait

