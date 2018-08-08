#!/bin/bash
IN=$1
THRES=${2:-0.01}
cat $IN \
| awk '/^[^#]/{ if($7<'$THRES') { print $5} }' \
| xargs -i iupac2meme {} 

#tail ame.tsv -n +2 | head  -n31 | cut -f5 | xargs -i iupac2meme {} >ame.meme
