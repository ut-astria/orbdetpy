#!/bin/bash

pushd $(dirname "$0") > /dev/null

export CLASSPATH=`find ../lib -name "*.jar" | tr "\n" ":"`

find . -name "*.java" | xargs javac

CLASS=`find . -name "*.class" | tr "\n" " "`
jar cvf ../lib/astria.jar $CLASS

popd > /dev/null
