/*
 *  (c) 2019 Open Source Geospatial Foundation - all rights reserved
 *  This code is licensed under the GPL 2.0 license, available at the root
 *  application directory.
 */

package org.geoserver.generatedgeometries;

import org.apache.wicket.markup.html.form.DropDownChoice;
import org.apache.wicket.model.IModel;
import org.apache.wicket.model.Model;
import org.geoserver.catalog.Catalog;
import org.geoserver.catalog.FeatureTypeInfo;
import org.geoserver.catalog.LayerInfo;
import org.geoserver.catalog.ResourcePool;
import org.geoserver.data.test.SystemTestData;
import org.geoserver.generatedgeometries.dummy.DummyGeometryGenerationStrategy;
import org.geoserver.generatedgeometries.longitudelatitude.LongLatTestData;
import org.geoserver.web.GeoServerApplication;
import org.geoserver.web.GeoServerWicketTestSupport;
import org.geoserver.web.data.resource.ResourceConfigurationPage;
import org.junit.Before;
import org.junit.Test;
import org.opengis.feature.type.FeatureType;

import java.io.IOException;
import java.util.Optional;

import static java.util.Collections.emptyMap;
import static org.geoserver.generatedgeometries.longitudelatitude.LongLatTestData.LONG_LAT_LAYER;
import static org.geoserver.generatedgeometries.longitudelatitude.LongLatTestData.LONG_LAT_NO_GEOM_ON_THE_FLY_LAYER;
import static org.geoserver.generatedgeometries.longitudelatitude.LongLatTestData.LONG_LAT_NO_GEOM_ON_THE_FLY_QNAME;
import static org.geoserver.generatedgeometries.longitudelatitude.LongLatTestData.LONG_LAT_QNAME;
import static org.geoserver.generatedgeometries.longitudelatitude.LongLatTestData.filenameOf;
import static org.geoserver.platform.TestGeoServerExtensionFinders.emptyFinder;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

public class GeneratedGeometryConfigurationPanelTest extends GeoServerWicketTestSupport {

    private static final String CREATE_GEOMETRY_LINK_PATH =
            "publishedinfo:tabs:panel:theList:0:content:content:createGeometryLink";

    private final GeoServerApplication application = mock(GeoServerApplication.class);

    @Before
    public void before() {
        System.setProperty("quietTests", "false");
    }

    @Override
    protected void onSetUp(SystemTestData testData) throws Exception {
        super.onSetUp(testData);
        setupLayerWithoutGeometry(testData);
        setupLayerWithGeometry(testData);
    }

    private void setupLayerWithoutGeometry(SystemTestData testData) throws IOException {
        testData.addVectorLayer(
                LONG_LAT_NO_GEOM_ON_THE_FLY_QNAME,
                emptyMap(),
                filenameOf(LONG_LAT_NO_GEOM_ON_THE_FLY_LAYER),
                LongLatTestData.class,
                getCatalog());
    }

    private void setupLayerWithGeometry(SystemTestData testData) throws IOException {
        Catalog catalog = getCatalog();
        testData.addVectorLayer(
                LONG_LAT_QNAME,
                emptyMap(),
                filenameOf(LONG_LAT_LAYER),
                LongLatTestData.class,
                catalog);
    }

    @Test
    public void testThatConfigurationIsNotAvailableForNonSimpleFeatureType() throws IOException {
        // given
        Catalog catalog = mock(Catalog.class);
        ResourcePool resourcePool = mock(ResourcePool.class);
        FeatureTypeInfo info = mock(FeatureTypeInfo.class);
        FeatureType featureType = mock(FeatureType.class);

        when(application.getCatalog()).thenReturn(catalog);
        when(catalog.getResourcePool()).thenReturn(resourcePool);
        when(resourcePool.getFeatureType(info)).thenReturn(featureType);

        IModel model = new Model<>(info);
        GeneratedGeometryConfigurationPanel panel =
                new GeneratedGeometryConfigurationPanel(
                        "id", model, emptyFinder(), () -> application);

        // when
        tester.startComponentInPage(panel);

        // then
        tester.assertVisible("id:content:errorMessage");
        IModel<?> defaultModel =
                tester.getComponentFromLastRenderedPage("id:content:errorMessage")
                        .getDefaultModel();
        assertEquals(defaultModel.getObject(), "incorrectFeatureType");
    }

    private DropDownChoice<GeometryGenerationStrategy> getStrategyDropDown() {
        return (DropDownChoice<GeometryGenerationStrategy>)
                tester.getComponentFromLastRenderedPage(
                        "publishedinfo:tabs:panel:theList:0:content:content:methodologyDropDown");
    }

    private Optional<? extends GeometryGenerationStrategy> findDummyStrategy(
            DropDownChoice<GeometryGenerationStrategy> dropDown) {
        return dropDown.getChoices()
                .stream()
                .filter(strategy -> strategy.getName().equals(DummyGeometryGenerationStrategy.NAME))
                .findFirst();
    }

    @Test
    public void testThatGeneratedGeometryPanelIsVisibleAndContainsStrategyWhenExtensionLoaded() {
        // given
        Catalog catalog = getGeoServerApplication().getCatalog();
        LayerInfo layerInfo = catalog.getLayerByName(getLayerId(LONG_LAT_NO_GEOM_ON_THE_FLY_QNAME));
        login();

        // when
        tester.startPage(new ResourceConfigurationPage(layerInfo, true));

        // then
        tester.assertComponent(
                "publishedinfo:tabs:panel:theList:0:content",
                GeneratedGeometryConfigurationPanel.class);
        DropDownChoice<GeometryGenerationStrategy> dropDown = getStrategyDropDown();
        Optional<? extends GeometryGenerationStrategy> dummyStrategyOptional =
                findDummyStrategy(dropDown);
        assertTrue(dummyStrategyOptional.isPresent());
    }

    @Test
    public void testThatCreateGeometryLinkActivatesSelectedStrategy() {
        // given
        Catalog catalog = getGeoServerApplication().getCatalog();
        LayerInfo layerInfo = catalog.getLayerByName(getLayerId(LONG_LAT_NO_GEOM_ON_THE_FLY_QNAME));
        login();
        tester.startPage(new ResourceConfigurationPage(layerInfo, true));
        DropDownChoice<GeometryGenerationStrategy> dropDown = getStrategyDropDown();
        GeometryGenerationStrategy dummyStrategy = findDummyStrategy(dropDown).get();
        dropDown.setModelObject(dummyStrategy);

        // when
        executeAjaxEventBehavior(CREATE_GEOMETRY_LINK_PATH, "click", "0");

        // then
        DummyGeometryGenerationStrategy strategy = (DummyGeometryGenerationStrategy) dummyStrategy;
        assertTrue(strategy.configured);
    }
}
