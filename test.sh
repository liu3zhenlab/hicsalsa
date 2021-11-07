check_sft () {
    pass=true
    cmd_status=`command -v $1`
    if [ -z $cmd_status ]; then
        pass=false
    fi  
    echo $pass
}


is_java=`check_sft java`
is_bwa=`check_sft xxx`
echo $is_java
echo $is_bwa

if $is_java && $is_bwa; then
	echo "yes"
fi
