#/bin/bash
main(){
	local IN=$1
	echo [IN]=$IN
	awk '{BNAME=FILENAME; sub("\\..*","",BNAME);   print $0 >> BNAME"."$8".bed"} END{print BNAME; }' $IN
}
main "$@"
