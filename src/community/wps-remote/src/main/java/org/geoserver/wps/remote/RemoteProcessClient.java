/* (c) 2014 Open Source Geospatial Foundation - all rights reserved
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */
package org.geoserver.wps.remote;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import javax.net.ssl.SSLContext;

import org.geoserver.catalog.Catalog;
import org.geoserver.catalog.DataStoreInfo;
import org.geoserver.catalog.LayerInfo;
import org.geoserver.catalog.StyleInfo;
import org.geoserver.catalog.WorkspaceInfo;
import org.geoserver.config.GeoServer;
import org.geoserver.importer.ImportContext;
import org.geoserver.importer.ImportTask;
import org.geoserver.importer.Importer;
import org.geoserver.importer.SpatialFile;
import org.geoserver.platform.GeoServerExtensions;
import org.geoserver.platform.GeoServerResourceLoader;
import org.opengis.feature.type.Name;
import org.opengis.util.ProgressListener;
import org.springframework.beans.factory.DisposableBean;

/**
 * Base class for the remote clients implementations. Those implementations will be plugged into GeoServer through the Spring app-context.
 * 
 * @author Alessio Fabiani, GeoSolutions
 * 
 */
public abstract class RemoteProcessClient implements DisposableBean {

    /** Whether this client is enabled or not from configuration */
    private boolean enabled;

    /** Whenever more instances of the client are available, they should be ordered by ascending priority */
    private int priority;

    /** The {@link RemoteProcessFactoryConfigurationWatcher} implementation */
    private final RemoteProcessFactoryConfigurationWatcher remoteProcessFactoryConfigurationWatcher;

    /** The registered {@link RemoteProcessFactoryListener} */
    private Set<RemoteProcessFactoryListener> remoteFactoryListeners = Collections
            .newSetFromMap(new ConcurrentHashMap<RemoteProcessFactoryListener, Boolean>());

    /** The registered {@link RemoteProcessClientListener} */
    private Set<RemoteProcessClientListener> remoteClientListeners = Collections
            .newSetFromMap(new ConcurrentHashMap<RemoteProcessClientListener, Boolean>());

    /**
     * The default Cosntructor
     * 
     * @param remoteProcessFactory
     */
    public RemoteProcessClient(
            RemoteProcessFactoryConfigurationWatcher remoteProcessFactoryConfigurationWatcher,
            boolean enabled, int priority) {
        this.remoteProcessFactoryConfigurationWatcher = remoteProcessFactoryConfigurationWatcher;
        this.enabled = enabled;
        this.priority = priority;
    }

    /**
     * @return the {@link RemoteProcessFactoryConfiguration} object
     */
    public RemoteProcessFactoryConfiguration getConfiguration() {
        return this.remoteProcessFactoryConfigurationWatcher.getConfiguration();
    }

    /**
     * Initialization method
     * 
     * @throws Exception
     */
    public abstract void init(SSLContext customSSLContext) throws Exception;

    /**
     * Destroy method
     * 
     * @throws Exception
     */
    public abstract void destroy() throws Exception;

    /**
     * @return the remoteFactoryListeners
     */
    public Set<RemoteProcessFactoryListener> getRemoteFactoryListeners() {
        return remoteFactoryListeners;
    }

    /**
     * @return the remoteClientListeners
     */
    public Set<RemoteProcessClientListener> getRemoteClientListeners() {
        return remoteClientListeners;
    }

    /**
     * @param enabled the enabled to set
     */
    public void setEnabled(boolean enabled) {
        this.enabled = enabled;
    }

    /**
     * Whether the plugin is enabled or not.
     * 
     * @return
     */
    public boolean isEnabled() {
        return this.enabled;
    }

    /**
     * @return the priority
     */
    public int getPriority() {
        return priority;
    }

    /**
     * @param priority the priority to set
     */
    public void setPriority(int priority) {
        this.priority = priority;
    }

    /**
     * Registers the {@link RemoteProcessFactoryListener} remoteClientListeners
     * 
     * @param listener
     */
    public void registerProcessFactoryListener(RemoteProcessFactoryListener listener) {
        remoteFactoryListeners.add(listener);
    }

    /**
     * De-registers the {@link RemoteProcessFactoryListener} remoteClientListeners
     * 
     * @param listener
     */
    public void deregisterProcessFactoryListener(RemoteProcessFactoryListener listener) {
        remoteFactoryListeners.remove(listener);
    }

    /**
     * Registers the {@link RemoteProcessClientListener} remoteClientListeners
     * 
     * @param listener
     */
    public void registerProcessClientListener(RemoteProcessClientListener listener) {
        remoteClientListeners.add(listener);
    }

    /**
     * De-registers the {@link RemoteProcessClientListener} remoteClientListeners
     * 
     * @param listener
     */
    public void deregisterProcessClientListener(RemoteProcessClientListener listener) {
        remoteClientListeners.remove(listener);
    }

    /**
     * Invoke the {@link RemoteProcessClient} execution
     * 
     * @param name
     * @param input
     * @param metadata
     * @param monitor
     * @return
     * @throws Exception
     */
    public abstract String execute(Name name, Map<String, Object> input,
            Map<String, Object> metadata, ProgressListener monitor) throws Exception;

    /**
     * Accessor for global geoserver instance from the test application context.
     */
    public GeoServer getGeoServer() {
        return (GeoServer) GeoServerExtensions.bean("geoServer");
    }

    /**
     * Accessor for global geoserver instance from the test application context.
     */
    public Importer getImporter() {
        return (Importer) GeoServerExtensions.bean("importer");
    }

    /**
     * 
     * @param wsName
     * @param dsName
     * @return
     */
    public DataStoreInfo createH2DataStore(String wsName, String dsName) {
        // create a datastore to import into
        Catalog cat = getGeoServer().getCatalog();

        WorkspaceInfo ws = wsName != null ? cat.getWorkspaceByName(wsName) : cat
                .getDefaultWorkspace();
        DataStoreInfo ds = cat.getFactory().createDataStore();
        ds.setWorkspace(ws);
        ds.setName(dsName);
        ds.setType("H2");

        GeoServerResourceLoader loader = cat.getResourceLoader();
        final String dataDir = loader.getBaseDirectory().getAbsolutePath();

        Map params = new HashMap();
        params.put("database", dataDir + "/" + dsName);
        params.put("dbtype", "h2");
        params.put("namespace", cat.getNamespaceByPrefix(ws.getName()).getURI());
        ds.getConnectionParameters().putAll(params);
        ds.setEnabled(true);
        cat.add(ds);

        return ds;
    }

    /**
     * @param value
     * @return
     * @throws IOException
     */
    public LayerInfo importLayer(File file, DataStoreInfo store, String defaultStyle,
            String targetWorkspace) throws Exception {
        Importer importer = getImporter();

        ImportContext context = (store != null ? importer.createContext(new SpatialFile(file),
                store) : importer.createContext(new SpatialFile(file)));

        if (context.getTasks() != null && context.getTasks().size() > 0) {
            ImportTask task = context.getTasks().get(0);

            if (targetWorkspace != null) {
                WorkspaceInfo ws = importer.getCatalog().getWorkspaceByName(targetWorkspace);
                if (ws != null) {
                    context.setTargetWorkspace(ws);
                } else {
                    context.setTargetWorkspace(importer.getCatalog().getDefaultWorkspace());
                }
            }
            
            if (defaultStyle != null) {
                StyleInfo style = importer.getCatalog().getStyleByName(defaultStyle);
                if (style == null && targetWorkspace != null) {
                    style = importer.getCatalog().getStyleByName(targetWorkspace, defaultStyle);
                }
                
                if (style != null) {
                    task.getLayer().setDefaultStyle(style);
                }
            }

            importer.run(context);

            if (context.getState() == ImportContext.State.COMPLETE) {
                if (context.getTasks() != null && context.getTasks().size() > 0) {
                    // ImportTask task = context.getTasks().get(0);
                    // assertEquals(ImportTask.State.READY, task.getState());

                    // assertEquals("the layer name", task.getLayer().getResource().getName());

                    task = context.getTasks().get(0);

                    return task.getLayer();
                }
            }
        }
        return null;
    }

}
