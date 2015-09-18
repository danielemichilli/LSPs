#!/bin/sh
HOST='lsps.co.nf'
USER='1948849_lsps'
PASSWD='J0140+56'

ftp -n $HOST <<END_SCRIPT
quote USER $USER
quote PASS $PASSWD
$1 $2

bye
END_SCRIPT
exit 0

