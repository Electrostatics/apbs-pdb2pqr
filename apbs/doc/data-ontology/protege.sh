#!/bin/sh

cd "/Applications/protege-5.0.0-beta-17-os-x/Protege-5.0.0-beta-17"

CMD_OPTIONS="-Dapple.laf.useScreenMenuBar=true -Dcom.apple.mrj.application.apple.menu.about.name=Protege -Xdock:name=Protege -Xdock:icon=Protege.icns"

jre/bin/java -Xmx500M -Xms250M \
     -Dlog4j.configuration=file:log4j.xml \
     -DentityExpansionLimit=100000000 \
     -Dfile.encoding=UTF-8 \
     -classpath bin/org.apache.felix.main.jar:bin/ProtegeLauncher.jar \
     $CMD_OPTIONS org.protege.osgi.framework.Launcher
#file:///Users/bake113/Git/ontologies/computing/ict.owl
