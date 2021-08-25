docker run --rm \
    -v $HOME:$HOME \
    -v /usr/local/bin/:/usr/local/bin \
    -v /tmp/spark:/tmp/spark \
    -v $HADOOP_HOME:$HADOOP_HOME \
    -w $HOME \
    broadinstitute/gatk:4.2.1.0 $HADOOP_HOME/bin/hadoop -h

HADOOP_HOME=/cluster/hadoop
SPARK_HOME=/cluster/spark
docker run --rm \
    -v $HOME:$HOME \
    -v /usr/local/helper_script/bin:/usr/local/helper_script/bin \
    -v /tmp/spark:/tmp/spark \
    -v $HADOOP_HOME/bin:$HADOOP_HOME/bin \
    --env HADOOP_HOME=$HADOOP_HOME \
    --env SPARK_HOME=$SPARK_HOME \
    -w $HOME \
    broadinstitute/gatk:4.2.1.0 /bin/bash -c "ls $HADOOP_HOME/bin/hadoop; ls $SPARK_HOME/bin/spark-submit; ls /usr/local/helper_script/bin"


docker run --rm \
    -v $HOME:$HOME \
    -v /usr/local/helper_script/bin:/usr/local/helper_script/bin \
    -v /tmp/spark:/tmp/spark \
    -v $HADOOP_HOME/bin:$HADOOP_HOME/bin \
    --env HADOOP_HOME=$HADOOP_HOME \
    --env SPARK_HOME=$SPARK_HOME \
    -w $HOME \
    broadinstitute/gatk:4.2.1.0 /bin/bash -c "$HADOOP_HOME/bin/hadoop -h; $SPARK_HOME/bin/spark-submit -h; ls /usr/local/helper_script/bin"

    docker run --rm \
    -v $HOME:$HOME \
    -v /usr/local/helper_script/bin:/usr/local/helper_script/bin \
    -v /tmp/spark:/tmp/spark \
    -v $HADOOP_HOME:$HADOOP_HOME \
    --env HADOOP_HOME=$HADOOP_HOME \
    --env SPARK_HOME=$SPARK_HOME \
    -w $HOME \
    broadinstitute/gatk:4.2.1.0 /usr/local/bin/yarn

    docker run --rm \
    -v /tmp/spark:/tmp/spark \
    -v /usr/local/bin:/usr/local/bin \
    broadinstitute/gatk:4.2.1.0 /usr/local/bin/yarn

    docker run -v $TOOLS:$TOOLS --env TOOLS=$TOOLS broadinstitute/gatk:4.2.1.0 ls $TOOLS