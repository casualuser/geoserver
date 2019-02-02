/*
 *  (c) 2019 Open Source Geospatial Foundation - all rights reserved
 *  This code is licensed under the GPL 2.0 license, available at the root
 *  application directory.
 */

package org.geoserver.generatedgeometries.longitudelatitude;

import org.geoserver.catalog.FeatureTypeInfo;
import org.geotools.data.FeatureSource;
import org.geotools.data.simple.SimpleFeatureSource;
import org.junit.Test;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.vfny.geoserver.global.ConfigurationException;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotSame;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;
import static org.mockito.BDDMockito.given;
import static org.mockito.Mockito.mock;

public class LongLatGeometryResourcePoolCallbackTest {

    private LongLatGeometryGenerationStrategy strategy =
            mock(LongLatGeometryGenerationStrategy.class);
    private LongLatGeometryResourcePoolCallback callback =
            new LongLatGeometryResourcePoolCallback(strategy);

    @Test
    public void testThatCanBuildFeatureType() {
        // given
        given(strategy.canHandle(null)).willReturn(false);
        FeatureTypeInfo info1 = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info1)).willReturn(true);
        FeatureTypeInfo info2 = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info2)).willReturn(false);

        // then
        assertFalse(callback.canBuildFeatureType(null, null));
        assertTrue(callback.canBuildFeatureType(info1, null));
        assertFalse(callback.canBuildFeatureType(info2, null));
    }

    @Test
    public void testThatReturnsSameFeatureTypeWhenCannotBuildIt() {
        // given
        FeatureTypeInfo info = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info)).willReturn(false);
        SimpleFeatureType featureType = mock(SimpleFeatureType.class);

        // when
        SimpleFeatureType builtFeatureType = callback.buildFeatureType(info, featureType);

        // then
        assertSame(featureType, builtFeatureType);
    }

    @Test
    public void testThatReturnsSameFeatureTypeInCaseOfErrorDuringBuilding()
            throws ConfigurationException {
        // given
        FeatureTypeInfo info = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info)).willReturn(true);
        SimpleFeatureType featureType = mock(SimpleFeatureType.class);
        given(strategy.defineGeometryAttributeFor(info, featureType))
                .willThrow(ConfigurationException.class);

        // when
        SimpleFeatureType builtFeatureType = callback.buildFeatureType(info, featureType);

        // then
        assertSame(featureType, builtFeatureType);
    }

    @Test
    public void testThatBuildsFeatureType() throws ConfigurationException {
        // given
        FeatureTypeInfo info = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info)).willReturn(true);
        SimpleFeatureType srcFeatureType = mock(SimpleFeatureType.class);
        SimpleFeatureType expectedFeatureType = mock(SimpleFeatureType.class);
        given(strategy.defineGeometryAttributeFor(info, srcFeatureType))
                .willReturn(expectedFeatureType);

        // when
        SimpleFeatureType builtFeatureType = callback.buildFeatureType(info, srcFeatureType);

        // then
        assertSame(expectedFeatureType, builtFeatureType);
    }

    @Test
    public void testThatCanWrapFeatureSource() {
        // given
        given(strategy.canHandle(null)).willReturn(false);
        FeatureTypeInfo info1 = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info1)).willReturn(true);
        FeatureTypeInfo info2 = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info2)).willReturn(false);

        // then
        assertFalse(callback.canWrapFeatureSource(null, null));
        assertTrue(callback.canWrapFeatureSource(info1, null));
        assertFalse(callback.canWrapFeatureSource(info2, null));
    }

    @Test
    public void testThatReturnsSameFeatureSourceWhenCannotWrapIt() {
        // given
        FeatureTypeInfo info = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info)).willReturn(false);
        FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = mock(FeatureSource.class);

        // when
        FeatureSource<SimpleFeatureType, SimpleFeature> wrapped =
                callback.wrapFeatureSource(info, featureSource);

        // then
        assertSame(featureSource, wrapped);
    }

    @Test
    public void testThatWrapsFeatureSource() {
        // given
        FeatureTypeInfo info = mock(FeatureTypeInfo.class);
        given(strategy.canHandle(info)).willReturn(true);
        FeatureSource<SimpleFeatureType, SimpleFeature> featureSource =
                mock(SimpleFeatureSource.class);

        // when
        FeatureSource<SimpleFeatureType, SimpleFeature> wrapped =
                callback.wrapFeatureSource(info, featureSource);

        // then
        assertNotSame(featureSource, wrapped);
        assertTrue(wrapped instanceof LongLatFeatureSource);
    }
}
