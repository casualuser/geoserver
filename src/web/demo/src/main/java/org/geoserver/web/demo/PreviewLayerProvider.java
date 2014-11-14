/* (c) 2014 Open Source Geospatial Foundation - all rights reserved
 * (c) 2001 - 2013 OpenPlans
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */
package org.geoserver.web.demo;


import java.util.ArrayList;
import java.util.Arrays;

import static org.geoserver.catalog.Predicates.*;

import java.util.Iterator;
import java.util.List;

import org.apache.wicket.extensions.markup.html.repeater.util.SortParam;
import org.apache.wicket.model.IModel;
import org.geoserver.catalog.Catalog;
import org.geoserver.catalog.LayerGroupInfo;
import org.geoserver.catalog.LayerInfo;
import org.geoserver.catalog.Predicates;
import org.geoserver.catalog.PublishedInfo;
import org.geoserver.catalog.util.CloseableIterator;
import org.geoserver.catalog.util.CloseableIteratorAdapter;
import org.geoserver.web.wicket.GeoServerDataProvider;
import org.opengis.filter.Filter;
import org.opengis.filter.sort.SortBy;

import com.google.common.base.Function;
import com.google.common.collect.Lists;


/**
 * Provides a filtered, sorted view over the catalog layers.
 * 
 * @author Andrea Aime - OpenGeo
 */
@SuppressWarnings("serial")
public class PreviewLayerProvider extends GeoServerDataProvider<PreviewLayer> {
    public static final Property<PreviewLayer> TYPE = new BeanProperty<PreviewLayer>(
            "type", "type");

    public static final AbstractProperty<PreviewLayer> NAME = new AbstractProperty<PreviewLayer>("name") {
        @Override
        public Object getPropertyValue(PreviewLayer item) {
            if (item.layerInfo != null) {
                return item.layerInfo.prefixedName();
            }
            if (item.groupInfo != null) {
                return item.groupInfo.prefixedName();
            }
            return null;
        }
    };

    public static final Property<PreviewLayer> TITLE = new BeanProperty<PreviewLayer>(
            "title", "title");
    
    public static final Property<PreviewLayer> ABSTRACT = new BeanProperty<PreviewLayer>(
            "abstract", "abstract", false);
    
    public static final Property<PreviewLayer> KEYWORDS = new BeanProperty<PreviewLayer>(
            "keywords", "keywords", false);

    public static final Property<PreviewLayer> COMMON = new PropertyPlaceholder<PreviewLayer>(
            "commonFormats");

    public static final Property<PreviewLayer> ALL = new PropertyPlaceholder<PreviewLayer>(
            "allFormats");

    public static final List<Property<PreviewLayer>> PROPERTIES = Arrays.asList(TYPE,
            NAME, TITLE, ABSTRACT, KEYWORDS, COMMON, ALL);
    
    @Override
    protected List<PreviewLayer> getItems() {
        // forced to implement this method as its abstract in the super class
        throw new UnsupportedOperationException(
                "This method should not be being called! "
                        + "We use the catalog streaming API");
    }

    @Override
    protected List<Property<PreviewLayer>> getProperties() {
        return PROPERTIES;
    }

    public IModel newModel(Object object) {
        return new PreviewLayerModel((PreviewLayer) object);
    }
    
    @Override
    public int size() {
        Filter filter = getFilter();
        return getCatalog().count(PublishedInfo.class, filter);
    }

    @Override
    public int fullSize() {
        Filter filter = Predicates.acceptAll();
        return getCatalog().count(PublishedInfo.class, filter);
    }
    
    @Override
    public Iterator<PreviewLayer> iterator(final int first, final int count) {
        Iterator<PreviewLayer> iterator = filteredItems(first, count);
        if (iterator instanceof CloseableIterator) {
            // don't know how to force wicket to close the iterator, lets return
            // a copy. Shouldn't be much overhead as we're paging
            try {
                return Lists.newArrayList(iterator).iterator();
            } finally {
                CloseableIteratorAdapter.close(iterator);
            }
        } else {
            return iterator;
        }
    }

    /**
     * Returns the requested page of layer objects after applying any keyword
     * filtering set on the page
     */
    @SuppressWarnings("resource")
    private Iterator<PreviewLayer> filteredItems(Integer first, Integer count) {
        final Catalog catalog = getCatalog();

        // global sorting
        final SortParam sort = getSort();
        final Property<PreviewLayer> property = getProperty(sort);

        SortBy sortOrder = null;
        if (sort != null) {
            if (property instanceof BeanProperty) {
                final String sortProperty = ((BeanProperty<PreviewLayer>) property)
                        .getPropertyPath();
                sortOrder = sortBy(sortProperty, sort.isAscending());
            } else if (property == NAME) {
                sortOrder = sortBy("prefixedName", sort.isAscending());
            }
        }

        Filter filter = getFilter();
        CloseableIterator<PublishedInfo> pi = catalog.list(PublishedInfo.class, filter, first,
                count, sortOrder);

        return CloseableIteratorAdapter.transform(pi, new Function<PublishedInfo, PreviewLayer>() {

            @Override
            public PreviewLayer apply(PublishedInfo input) {
                if (input instanceof LayerInfo) {
                    return new PreviewLayer((LayerInfo) input);
                } else if (input instanceof LayerGroupInfo) {
                    return new PreviewLayer((LayerGroupInfo) input);
                }
                return null;
            }
        });
    }
    
    @Override
    protected Filter getFilter() {
        Filter filter = super.getFilter();

        // need to get only advertised and enabled layers
        Filter enabledFilter = Predicates.or(Predicates.isNull("resource.enabled"),
                Predicates.equal("resource.enabled", true));
        Filter advertisedFilter = Predicates.or(Predicates.isNull("resource.advertised"),
                Predicates.equal("resource.advertised", true));

        return Predicates.and(filter, enabledFilter, advertisedFilter);
    }


}
