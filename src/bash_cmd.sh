#!/bin/sh
HOST='lsps.co.nf'
USER='1948849_lsps'
PASSWD='J0140+56'

ftp -n $HOST <<END_SCRIPT
quote USER $USER
quote PASS $PASSWD

if [ "$1" = "cd" ]; then
  ${@:1:2}
  ${@:3}
else
  $@
fi

bye
END_SCRIPT
exit 0

