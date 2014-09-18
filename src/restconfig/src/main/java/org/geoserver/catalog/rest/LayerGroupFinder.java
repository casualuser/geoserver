/* (c) 2014 Open Source Geospatial Foundation - all rights reserved
 * (c) 2001 - 2013 OpenPlans
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */
package org.geoserver.catalog.rest;

import org.geoserver.catalog.Catalog;
import org.geoserver.rest.RestletException;
import org.restlet.data.Form;
import org.restlet.data.Method;
import org.restlet.data.Request;
import org.restlet.data.Response;
import org.restlet.data.Status;
import org.restlet.resource.Resource;

public class LayerGroupFinder extends AbstractCatalogFinder {

    public LayerGroupFinder(Catalog catalog) {
        super(catalog);
    }
    
    @Override
    public Resource findTarget(Request request, Response response) {
        String ws = getAttribute(request, "workspace");
        String lg = getAttribute(request, "layergroup");
        
        if (ws != null && catalog.getWorkspaceByName(ws) == null) {
            throw new RestletException("No such workspace: " + ws, Status.CLIENT_ERROR_NOT_FOUND);
        }

        if ( lg != null) { 
            Form form = request.getResourceRef().getQueryAsForm();
            String quietOnNotFoundS=form.getFirstValue("quietOnNotFound", "False");
            boolean quietOnNotFound=Boolean.parseBoolean(quietOnNotFoundS);            
            //ensure it exists
            if (ws != null && catalog.getLayerGroupByName( ws, lg ) == null) {
                if(quietOnNotFound){
                    return null;
                }
                throw new RestletException(String.format("No such layer group %s in workspace %s", 
                    lg, ws), Status.CLIENT_ERROR_NOT_FOUND );
            }
            if (ws == null && catalog.getLayerGroupByName( lg ) == null) {
                if(quietOnNotFound){
                    return null;
                }
                throw new RestletException( "No such layer group " + lg, Status.CLIENT_ERROR_NOT_FOUND );
            }
        }

        if ( lg == null && request.getMethod() == Method.GET ) {
            return new LayerGroupListResource( getContext(), request, response, catalog );
        }
        
        return new LayerGroupResource( getContext(), request, response, catalog );
    }

}
