/* (c) 2014 Open Source Geospatial Foundation - all rights reserved
 * (c) 2001 - 2013 OpenPlans
 * This code is licensed under the GPL 2.0 license, available at the root
 * application directory.
 */
package org.geoserver.wps.remote.plugin;

import java.io.UnsupportedEncodingException;
import java.net.URLDecoder;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

import org.geoserver.wps.remote.RemoteProcessClientListener;
import org.geotools.util.logging.Logging;
import org.jivesoftware.smack.packet.Message;
import org.jivesoftware.smack.packet.Packet;

/**
 * 
 * 
 * @author Alessio Fabiani, GeoSolutions
 * 
 */
public class XMPPErrorMessage implements XMPPMessage {

    /** The LOGGER */
    public static final Logger LOGGER = Logging.getLogger(XMPPMessage.class.getPackage().getName());

    @Override
    public boolean canHandle(Map<String, String> signalArgs) {
        if (signalArgs != null && signalArgs.get("topic") != null)
            return signalArgs.get("topic").equals("error");
        return false;
    }

    @Override
    public void handleSignal(XMPPClient xmppClient, Packet packet, Message message,
            Map<String, String> signalArgs) {

        Map<String, Object> metadata = new HashMap<String, Object>();
        metadata.put("serviceJID", packet.getFrom());

        Exception cause = null;
        try {
            cause = new Exception(URLDecoder.decode(signalArgs.get("message"), "UTF-8"));
        } catch (UnsupportedEncodingException e) {
            cause = e;
        }
        final String pID = (signalArgs != null ? signalArgs.get("id") : null);

        // NOTIFY SERVICE
        final String serviceJID = message.getFrom();
        xmppClient.sendMessage(serviceJID, "topic=abort");

        // NOTIFY LISTENERS
        for (RemoteProcessClientListener listener : xmppClient.getRemoteClientListeners()) {
            listener.exceptionOccurred(pID, cause, metadata);
        }

    }

}
