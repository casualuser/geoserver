package org.geoserver.wps.jdbc;
/* (c) 2015 Open Source Geospatial Foundation - all rights reserved
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.StringReader;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.parsers.ParserConfigurationException;

import org.geoserver.wps.ProcessStatusStore;
import org.geoserver.wps.executor.ExecutionStatus;
import org.geoserver.wps.executor.ProcessState;
import org.geoserver.wps.xml.WPSConfiguration;
import org.geotools.data.DataStore;
import org.geotools.data.DataUtilities;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.Query;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.NameImpl;
import org.geotools.feature.SchemaException;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.filter.text.cql2.CQLException;
import org.geotools.filter.text.ecql.ECQL;
import org.geotools.util.logging.Logging;
import org.geotools.wps.WPS;
import org.geotools.xml.Encoder;
import org.geotools.xml.Parser;
import org.opengis.feature.Property;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.feature.type.Name;
import org.opengis.filter.Filter;
import org.xml.sax.SAXException;

import net.opengis.ows11.Ows11Factory;
import net.opengis.wps10.ExecuteType;
import net.opengis.wps10.Wps10Factory;

/**
 * A class that stores WPS session statuses in JDBC datastores.
 * 
 * @author Ian Turton
 */

public class JDBCStatusStore implements ProcessStatusStore {

    private static final String STACKTRACESEPERATOR = "/";

    static final Logger LOGGER = Logging.getLogger(JDBCStatusStore.class);

    static final Wps10Factory WPSFACTORY = Wps10Factory.eINSTANCE;

    static final Ows11Factory OWSFACTORY = Ows11Factory.eINSTANCE;

    static final String STACK_TRACE = "stackTrace";

    static final String EXCEPTION_MESSAGE = "exceptionMessage";

    static final String EXCEPTION_CLASS = "exceptionClass";

    static final String ASYNC = "async";

    static final String USER_NAME = "user";

    static final String TASK = "task";

    static final String SIMPLE_PROCESS_NAME = "processName";

    static final String PROPERTIES = "properties";

    static final String PROGRESS = "progress";

    static final String PROCESS_NAME_URI = "processNameURI";

    static final String PROCESS_NAME = "processNameImpl";

    static final String PHASE = "phase";

    static final String NODE_ID = "node";

    static final String COMPLETION = "completionTime";

    static final String LASTUPDATE = "lastUpdated";

    static final String CREATION = "creationTime";

    static final String PROCESS_ID = "processId";

    static final String STATUS = "status";

    static final String EXECUTION_ID = "exceptionId";

    private static final String REQUEST = "request";

    DataStore statuses;

    static SimpleFeatureType schema;

    public JDBCStatusStore(JDBCStatusStoreLoader loader) {
        this(loader.getStore());
    }

    public JDBCStatusStore(DataStore store) {
        if (store == null) {
            throw new RuntimeException(
                    "Attempted to create a JDBCStatusStore with a null datastore");
        }
        statuses = store;
        try {
            boolean statusSchema = statuses.getNames().contains(new NameImpl(STATUS));
            if (statusSchema) {
                LOGGER.fine("using exisiting status store DB " + STATUS);
                SimpleFeatureSource source = statuses.getFeatureSource(STATUS);
                schema = source.getSchema();
            } else {
                LOGGER.fine("creating new DB table for statuses");
                try {
                  //@formatter:off
                    schema = DataUtilities.createType(STATUS,
                              PROCESS_ID+":String," 
                            + CREATION+":Date," 
                            + LASTUPDATE+":Date,"
                            + COMPLETION+":Date," 
                            + NODE_ID+":String," 
                            + PHASE+":String,"
                            + PROCESS_NAME+":String," 
                            + PROCESS_NAME_URI+":String,"
                            + PROGRESS+":Float," 
                            + REQUEST+":String," 
                            + PROPERTIES+":String," 
                            + SIMPLE_PROCESS_NAME+":String,"
                            + TASK+":String," 
                            + USER_NAME+":String," 
                            + ASYNC+":Boolean,"
                            + EXCEPTION_CLASS+":String,"
                            + EXCEPTION_MESSAGE+":String," 
                            + STACK_TRACE+":String");
                } catch (SchemaException e) {
                   LOGGER.fine(e.getLocalizedMessage());
                }
              //@formatter:on
                statuses.createSchema(schema);
            }
        } catch (IOException e1) {
            LOGGER.info("failed to create WPS status store");
            LOGGER.fine(e1.getMessage());
        }

    }

    @Override
    public void save(ExecutionStatus status) {
        DefaultTransaction transaction = new DefaultTransaction("create");
        try {
            SimpleFeatureSource source = statuses.getFeatureSource(STATUS);
            if (source instanceof SimpleFeatureStore) {
                SimpleFeatureStore store = (SimpleFeatureStore) source;
                store.setTransaction(transaction);
                SimpleFeature feature = statusToFeature(status);
                FeatureCollection<SimpleFeatureType, SimpleFeature> featureCollection = DataUtilities
                        .collection(feature);
                // if the feature exists delete it
                Filter filter = ECQL.toFilter(PROCESS_ID + " = '" + status.getExecutionId() + "'");
                store.removeFeatures(filter);
                store.addFeatures(featureCollection);
                transaction.commit();
            }
        } catch (IOException e) {
            LOGGER.info("failed to save WPS Status " + status.getExecutionId());
            LOGGER.fine(e.getMessage());

            try {
                transaction.rollback();
            } catch (IOException e1) {
                // everything has gone wrong?
                LOGGER.fine("failed to rollback failed save");
            }
        } catch (CQLException e) {
            // should really never be thrown
        } finally {
            transaction.close();
        }

    }

    @Override
    public ExecutionStatus get(String executionId) {
        LOGGER.fine("getting status " + executionId);
        try {
            SimpleFeatureSource source = statuses.getFeatureSource(STATUS);

            Filter filter = ECQL.toFilter(PROCESS_ID + " = '" + executionId + "'");
            SimpleFeatureCollection features = source.getFeatures(filter);
            SimpleFeature f = DataUtilities.first(features);
            ExecutionStatus stat = featureToStatus(f);
            return stat;
        } catch (IOException | CQLException e) {
            LOGGER.info("failed to fetch WPS Status " + executionId);
            LOGGER.fine(e.getMessage());
        }
        return null;
    }

    @Override
    public ExecutionStatus remove(String executionId) {
        LOGGER.fine("removing status " + executionId);
        try {
            SimpleFeatureSource source = statuses.getFeatureSource(STATUS);

            Filter filter = ECQL.toFilter(PROCESS_ID + " = '" + executionId + "'");
            SimpleFeatureCollection features = source.getFeatures(filter);
            SimpleFeature f = DataUtilities.first(features);
            ExecutionStatus stat = featureToStatus(f);
            DefaultTransaction transaction = new DefaultTransaction("create");
            try {
                if (source instanceof SimpleFeatureStore) {
                    SimpleFeatureStore store = (SimpleFeatureStore) source;
                    store.setTransaction(transaction);
                    store.removeFeatures(filter);
                    transaction.commit();
                }
            } catch (IOException e) {
                LOGGER.info("failed to remove WPS Status " + executionId);
                LOGGER.fine(e.getMessage());
                try {
                    transaction.rollback();
                } catch (IOException e1) {
                    // don't care what happens here
                }
            } finally {
                transaction.close();
            }
            return stat;
        } catch (IOException | CQLException e) {
            LOGGER.info("failed to remove WPS Status " + executionId);
            LOGGER.fine(e.getMessage());

        }
        return null;
    }

    @Override
    public int remove(Filter filter) {
        LOGGER.fine("removing statuses matching " + filter);
        int ret = 0;
        DefaultTransaction transaction = new DefaultTransaction("create");
        try {
            SimpleFeatureSource source = statuses.getFeatureSource(STATUS);
            SimpleFeatureCollection features = source.getFeatures(filter);
            ret = features.size();
            if (ret == 0) {
                return ret;
            }
            SimpleFeatureIterator itr = features.features();

            itr.close();
            if (source instanceof SimpleFeatureStore) {
                SimpleFeatureStore store = (SimpleFeatureStore) source;
                store.setTransaction(transaction);
                store.removeFeatures(filter);
                transaction.commit();
            } else {
                LOGGER.info("Readonly status store found, probably not what you wanted");
            }
        } catch (IOException e) {
            LOGGER.info("failed to remove WPS Status matching " + filter);
            LOGGER.fine(e.getMessage());
            try {
                transaction.rollback();
            } catch (IOException e1) {
                // don't care what happens here
            }
        } finally {
            transaction.close();
        }

        return ret;
    }

    @Override
    public List<ExecutionStatus> list(Query query) {
        LOGGER.fine("listing statuses matching " + query);
        SimpleFeatureSource source;
        try {
            ArrayList<ExecutionStatus> ret = new ArrayList<>();
            source = statuses.getFeatureSource(STATUS);
            LOGGER.fine("requesting " + query);

            SimpleFeatureCollection features = source.getFeatures(query);
            int size = features.size();
            LOGGER.fine(query + " " + size);

            SimpleFeatureIterator itr = features.features();
            try {
                while (itr.hasNext()) {
                    SimpleFeature f = itr.next();
                    if (LOGGER.isLoggable(Level.FINE)) {
                        LOGGER.fine("adding " + f);
                    }
                    ret.add(featureToStatus(f));
                }
                return ret;
            } finally {
                itr.close();
            }
        } catch (IOException e) {
            LOGGER.info("failed to fetch any WPS statuses");
            LOGGER.fine(e.getMessage());
        }
        return Collections.emptyList();
    }

    protected SimpleFeature statusToFeature(ExecutionStatus status) {
        SimpleFeatureBuilder builder = new SimpleFeatureBuilder(schema);
        builder.set(PROCESS_ID, status.getExecutionId());
        builder.set(CREATION, status.getCreationTime());
        builder.set(LASTUPDATE, status.getLastUpdated());
        builder.set(COMPLETION, status.getCompletionTime());
        builder.set(NODE_ID, status.getNodeId());
        builder.set(PHASE, status.getPhase());
        Name processName = status.getProcessName();
        builder.set(PROCESS_NAME, processName.getLocalPart());
        builder.set(PROCESS_NAME_URI, processName.getNamespaceURI());

        builder.set(PROGRESS, status.getProgress());
        ExecuteType request = status.getRequest();
        if (request != null) {
            builder.set(REQUEST, serializeRequest(request));
        }
        builder.set(SIMPLE_PROCESS_NAME, status.getSimpleProcessName());
        builder.set(TASK, status.getTask());
        builder.set(USER_NAME, status.getUserName());
        builder.set(ASYNC, status.isAsynchronous());
        Throwable exception = status.getException();
        if (exception != null) {
            builder.set(EXCEPTION_CLASS, exception.getClass().getName());
            builder.set(EXCEPTION_MESSAGE, exception.getMessage());
            StackTraceElement[] stackTrace = exception.getStackTrace();
            StringBuffer buf = new StringBuffer();
            for (StackTraceElement el : stackTrace) {
                buf.append(el.getClassName()).append(STACKTRACESEPERATOR);
                buf.append(el.getFileName()).append(STACKTRACESEPERATOR);
                buf.append(el.getMethodName()).append(STACKTRACESEPERATOR);
                buf.append(el.getLineNumber());
                buf.append("\n");
            }

            builder.set(STACK_TRACE, buf.toString());
        }
        SimpleFeature feature = builder.buildFeature(null);
        return feature;
    }

    private String serializeRequest(ExecuteType request) {
        String ret = "";
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        Encoder e = new Encoder(new WPSConfiguration());
        e.setIndenting(true);
        try {
            e.encode(request, WPS.Execute, out);
            ret = out.toString();
        } catch (IOException e1) {
            LOGGER.info("Problem encountered encoding WPS Request");
        }

        return ret;

    }

    protected ExecutionStatus featureToStatus(SimpleFeature f) {
        HashMap<String, Object> attrs = new HashMap<>();
        if (f == null) {
            return null;
        }
        for (Property p : f.getProperties()) {
            //
            if (p.getValue() != null) {
                attrs.put(p.getName().toString(), p.getValue());
            }
        }
        Name processName = new NameImpl((String) attrs.get(PROCESS_NAME));
        String executionId = (String) attrs.get(PROCESS_ID);
        boolean asynchronous = (boolean) attrs.get(ASYNC);

        ExecutionStatus status = new ExecutionStatus(processName, executionId, asynchronous);
        if (attrs.containsKey(REQUEST)) {
            ExecuteType request = buildRequest(attrs);

            status.setRequest(request);
        }
        String phase = (String) attrs.get(PHASE);
        ProcessState state = ProcessState.valueOf(phase);
        status.setPhase(state);

        status.setProgress((float) attrs.get(PROGRESS));
        status.setTask((String) attrs.get(TASK));
        status.setUserName((String) attrs.get(USER_NAME));
        if (attrs.containsKey(EXCEPTION_MESSAGE)) {
            status.setException(buildException(attrs, status));
        }
        // set the time stamps last as other items may set them for us but
        // we know best here!
        status.setCreationTime((Date) attrs.get(CREATION));
        if (attrs.containsKey(COMPLETION)) {
            status.setCompletionTime((Date) attrs.get(COMPLETION));
        }
        if (attrs.containsKey(LASTUPDATE)) {
            status.setLastUpdated((Date) attrs.get(LASTUPDATE));
        }

        return status;
    }

    private ExecuteType buildRequest(HashMap<String, Object> attrs) {
        String req = (String) attrs.get(REQUEST);
        BufferedReader reader = new BufferedReader(new StringReader(req));
        org.geotools.xml.Parser parser = new Parser(new WPSConfiguration());
        ExecuteType request = null;
        try {
            request = (ExecuteType) parser.parse(reader);
        } catch (IOException | SAXException | ParserConfigurationException e) {
            LOGGER.info("Problem building WPS request for status");
            LOGGER.fine(e.getMessage());
        }

        return request;
    }

    private Exception buildException(HashMap<String, Object> attrs, ExecutionStatus status) {
        String message = (String) attrs.get(EXCEPTION_MESSAGE);
        Exception exc = new Exception(message);
        // see if we can rebuild the exception
        try {
            Constructor<?> con = this.getClass().getClassLoader()
                    .loadClass((String) attrs.get(EXCEPTION_CLASS)).getConstructor(String.class);
            exc = (Exception) con.newInstance(message);

        } catch (InstantiationException | IllegalAccessException | ClassNotFoundException
                | NoSuchMethodException | SecurityException | IllegalArgumentException
                | InvocationTargetException e) {
            // too bad, I don't care
            LOGGER.log(Level.FINE,"COuldn't reinstaniate Exception for WPS status",e);
        }

        String r = (String) attrs.get(STACK_TRACE);
        ArrayList<StackTraceElement> trace = new ArrayList<>();

        for (String line : r.split("\n")) {
            String[] parts = line.split(STACKTRACESEPERATOR);
            String declaringClass = parts[0];
            String fileName = parts[1];
            String methodName = parts[2];
            int lineNumber = Integer.parseInt(parts[3]);
            StackTraceElement t = new StackTraceElement(declaringClass, methodName, fileName,
                    lineNumber);

            trace.add(t);
        }
        exc.setStackTrace(trace.toArray(new StackTraceElement[] {}));
        return exc;
    }

    @Override
    public boolean supportsPredicate() {
        return false;
    }

    @Override
    public boolean supportsPaging() {
        return true;
    }

}
