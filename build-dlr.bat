@echo off
mvn -f src/pom.xml clean install -Pwps,wps-cluster-hazelcast,wps-remote,importer,security,dyndimension,colormap,netcdf,netcdf-out,rest-ext,jms-cluster -DskipTests