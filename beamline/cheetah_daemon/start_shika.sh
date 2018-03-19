#!/bin/sh
export PHENIX_TRUST_OTHER_ENV=1
yamtbx.shika_backend nproc=0 bl=32xu > ~/shika_backend.log 2>&1 &
/oys/xtal/cheetah/eiger-zmq/bin/cheetah.streaming_client \
 --bl=32xu --nproc=64 --logdir=/isilon/cluster/log/shika \
 --eiger-host=192.168.215.43:9999 \
 --vent-hosts="" \
 --result-host="127.0.0.1:5558" \
 --pub-host="" \
 > ~/shika_cheetah.log 2>&1 &
wait

