/* (c) 2016 Open Source Geospatial Foundation - all rights reserved
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */

package org.geoserver.netcdf;

import static org.junit.Assert.*;

import org.junit.Test;

public class NetCDFStatusTest {

    @Test
    public void testNetCDFStatus() {
        NetCDFStatus status = new NetCDFStatus();

        // these should always return "something"
        assertTrue(status.getModule().length() > 0);
        assertTrue(status.getName().length() > 0);
        assertTrue(status.getComponent().get().length() > 0);
        assertTrue(status.getMessage().get().length() > 0);

        assertTrue(status.getVersion().isPresent());
        assertTrue(status.getMessage().get().contains("NETCDF-4 Binary"));
    }
}
