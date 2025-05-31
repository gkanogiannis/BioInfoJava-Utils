/*
 *
 * BioInfoJava-Utils 
 *
 * Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */
package com.gkano.bioinfo.javautils;

import java.io.IOException;
import java.util.Properties;

public class VersionInfo {
    public static final String VERSION;
    public static final String BUILD_TIME;

    static {
        Properties props = new Properties();
        try {
            props.load(VersionInfo.class.getClassLoader().getResourceAsStream("version.properties"));
        } catch (IOException | NullPointerException e) {
            System.err.println("Warning: version.properties not found or unreadable.");
        }
        VERSION = props.getProperty("version", "unknown");
        BUILD_TIME = props.getProperty("buildTime", "unknown");
    }
}