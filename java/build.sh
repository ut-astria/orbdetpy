#!/bin/bash

pushd $(dirname "$0") > /dev/null

export CLASSPATH=`find ../lib -name "*.jar" | tr "\n" ":"`

find . -name "*.class" | xargs rm
find . -name "*.java" | xargs javac -Xlint:unchecked

CLASS=`find . -name "*.class" | tr "\n" " "`
jar cf ../lib/astria.jar $CLASS

popd > /dev/null
