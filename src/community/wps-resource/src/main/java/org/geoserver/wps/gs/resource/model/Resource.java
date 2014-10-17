/* (c) 2014 Open Source Geospatial Foundation - all rights reserved
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */
package org.geoserver.wps.gs.resource.model;

import java.util.logging.Logger;

import org.geoserver.wps.gs.resource.model.translate.TranslateContext;
import org.geotools.util.logging.Logging;

/**
 * This is the base class for each WPS Resource concrete implementation.
 * 
 * @author alessio.fabiani
 * 
 */
public abstract class Resource {

    static protected Logger LOGGER = Logging.getLogger(Resource.class);

    private String type;

    protected String name;

    protected boolean persistent;

    protected TranslateContext translateContext;

    /**
     * @return the type
     */
    public String getType() {
        return type;
    }

    /**
     * @param type the type to set
     */
    public void setType(String type) {
        this.type = type;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }

    /**
     * @param name the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the persistent
     */
    public boolean isPersistent() {
        return persistent;
    }

    /**
     * @param persistent the persistent to set
     */
    public void setPersistent(boolean persistent) {
        this.persistent = persistent;
    }

    /**
     * @return the translateContext
     */
    public TranslateContext getTranslateContext() {
        return translateContext;
    }

    /**
     * @param translateContext the translateContext to set
     */
    public void setTranslateContext(TranslateContext context) {
        this.translateContext = context;
    }

    /**
     * 
     * @return
     */
    public boolean isWellDefined() {
        boolean res = true;

        if (getType() == null || getName() == null)
            return false;

        return res && resourcePropertiesConsistencyCheck();
    }

    /**
     * 
     * @return
     */
    protected abstract boolean resourcePropertiesConsistencyCheck();
}
