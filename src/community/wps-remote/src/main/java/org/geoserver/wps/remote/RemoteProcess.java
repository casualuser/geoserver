/* (c) 2014 Open Source Geospatial Foundation - all rights reserved
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */
package org.geoserver.wps.remote;

import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.CountDownLatch;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.geotools.process.Process;
import org.geotools.process.ProcessException;
import org.geotools.util.logging.Logging;
import org.opengis.feature.type.Name;
import org.opengis.util.ProgressListener;

/**
 * Stub for the remote processes generated at run-time by the {@link RemoteProcessFactory} upon a {@link RemoteProcessClient} registration request.
 * 
 * @author Alessio Fabiani, GeoSolutions
 * 
 */
public class RemoteProcess implements Process, RemoteProcessClientListener {

    /** The LOGGER. */
    private static final Logger LOGGER = Logging.getLogger(RemoteProcess.class);

    /** The Process Name; declared by the remote service */
    private Name name;

    /** The {@link RemoteProcessClient} */
    private RemoteProcessClient remoteClient;

    /** A generic kvp map containing client specific implementation properties */
    private Map<String, Object> metadata;

    /** The Process Outputs; declared by the remote service */
    private Map<String, Object> outputs;

    /** Whether the Process is still running or not */
    private boolean running;

    /**
     * A Process ID generated by the {@link RemoteProcessClient}; this is used to uniquely identify the remote service sending commands and messages
     * to this {@link RemoteProcess} instance
     */
    private String pid;

    /** The progess listrener. */
    private ProgressListener listener;

    /**
     * Whether the remote service raised and exception or not. This property contains the cause and is instantiated by the {@link RemoteProcessClient}
     */
    private Exception exception;

    /** The semaphore */
    CountDownLatch doneSignal = new CountDownLatch(1);

    /**
     * Constructs a new stub for the {@link RemoteProcess} execution. Metadata is a kvp map containing specific properties of the
     * {@link RemoteProcessClient} instance
     * 
     * @param name
     * @param remoteClient
     * @param metadata
     */
    public RemoteProcess(Name name, RemoteProcessClient remoteClient,
            Map<String, Object> metadata) {
        this.name = name;
        this.remoteClient = remoteClient;
        this.metadata = metadata;
    }

    @Override
    public Map<String, Object> execute(Map<String, Object> input, ProgressListener monitor) {

        try {
            // Generating a unique Process ID

            LOGGER.info("Generating a unique Process ID for Remote Process [" + name
                    + "] with the following parameters:");
            LOGGER.info(" - name: " + name);
            LOGGER.info(" - input: " + input);
            LOGGER.info(" - metadata: " + metadata);
            LOGGER.info(" - monitor: " + monitor);

            if (remoteClient == null) {
                LOGGER.log(Level.SEVERE, "Cannot execute Remote Process [" + name
                        + "] since the RemoteClient is not available!");
                throw new Exception("Cannot execute Remote Process [" + name
                        + "] since the RemoteClient is not available!");
            }
            listener = monitor;
            pid = remoteClient.execute(name, input, metadata, monitor);
            LOGGER.info("Starting the execution of Remote Process with pId [" + pid + "]");
            running = pid != null;
            if (running && (listener != null && !listener.isCanceled())) {
                remoteClient.registerProcessClientListener(this);

                // doneSignal.await(timeout, unit); // TIMEOUT TODO
                doneSignal.await();
            }
            LOGGER.info("Stopping the execution of Remote Process with pId [" + pid + "]");
        } catch (Exception e) {
            if (listener != null) {
                listener.exceptionOccurred(e);
            }
            LOGGER.log(Level.SEVERE, "The Remote Process with pId [" + pid + "] rasied an Exeption",
                    e);
            throw new ProcessException(e);
        } finally {
            remoteClient.deregisterProcessClientListener(this);
        }

        // forward the Exception if necessary
        if (exception != null) {
            LOGGER.log(Level.SEVERE, "The Remote Service associated to the Process with pId [" + pid
                    + "] rasied an Exeption", exception);
            throw new ProcessException(exception);
        }

        // check if the Process has been cancelled
        if (listener != null && listener.isCanceled()) {
            LOGGER.log(Level.WARNING, "The Remote Service associated to the Process with pId ["
                    + pid + "] has been cancelled");
            throw new ProcessException("The Remote Service associated to the Process with pId ["
                    + pid + "] has been cancelled");
        }

        return outputs;
    }

    /**
     * @return the running
     */
    public boolean isRunning() {
        return running;
    }

    /**
     * @param running the running to set
     */
    public void setRunning(boolean running) {
        this.running = running;
    }

    /**
     * @return the outputs
     */
    public Map<String, Object> getOutputs() {
        return outputs;
    }

    /**
     * @param outputs the outputs to set
     */
    public void setOutputs(Map<String, Object> outputs) {
        this.outputs = outputs;
    }

    @Override
    public void progress(final String pId, final Double progress) {
        if (pId.equals(pid)) {
            listener.progress(progress.floatValue());
        }

        if (listener.isCanceled()) {
            doneSignal.countDown();
        }
    }

    @Override
    public void complete(String pId, Object outputs) {
        if (pId.equals(pid)) {
            listener.complete();

            try {
                this.outputs = (Map<String, Object>) outputs;
            } catch (Exception e) {
                exception = e;
                LOGGER.log(Level.SEVERE,
                        "The Remote Service associated to the Process with pId [" + pid
                                + "] rasied an Exeption while setting the outputs on completion",
                        exception);
                this.outputs = null;
            }

            running = false;
            doneSignal.countDown();
        }
    }

    @Override
    public void exceptionOccurred(final String pId, Exception cause, Map<String, Object> metadata) {
        if (pId != null && pId.equals(pid)) {
            listener.exceptionOccurred(cause);
            exception = cause;
            running = false;
        } else if (metadata != null) {
            boolean metadataIsEqual = true;

            for (Entry<String, Object> entry : metadata.entrySet()) {
                if (!this.metadata.containsKey(entry.getKey())
                        || this.metadata.get(entry.getKey()) != entry.getValue()) {
                    metadataIsEqual = false;
                    break;
                }
            }

            if (metadataIsEqual) {
                listener.exceptionOccurred(cause);
                exception = cause;
                running = false;
            }
        }
        doneSignal.countDown();
    }

}
