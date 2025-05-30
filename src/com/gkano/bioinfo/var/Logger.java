package com.gkano.bioinfo.var;

import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class Logger {

    private static final SimpleDateFormat DATE_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private static boolean verbose = false;
    private static boolean lastWasCarret = false;

    public static void setVerbose(boolean v) {
        verbose = v;
    }

    public static void info(Object caller, String msg) {
        if (verbose) {
            info(caller, msg, System.err);
        }
    }

    public static void info(Object caller, String msg, PrintStream ps) {
        if (lastWasCarret) {
            ps.println(); 
            lastWasCarret = false;
        }
        if (verbose) {
            ps.println("[INFO]\t" + caller.getClass().getSimpleName() + "\t" + timestamp() + " - " + msg);
            ps.flush();
        }
    }

    public static void infoCarret(Object caller, String msg) {
        if (verbose) {
            infoCarret(caller, msg, System.err);
        }
    }

    public static void infoCarret(Object caller, String msg, PrintStream ps) {
        if (verbose) {
            ps.print("\033[2K\r" +"[INFO]\t" + caller.getClass().getSimpleName() + "\t" + timestamp() + " - " + msg);
            ps.flush();
            lastWasCarret = true;
        }
    }

    public static void warn(Object caller, String msg) {
        warn(caller, msg, System.err);
    }

    public static void warn(Object caller, String msg, PrintStream ps) {
        ps.println("[WARN]\t" + caller.getClass().getSimpleName() + "\t" + timestamp() + " - " + msg);
        ps.flush();
    }

    public static void error(Object caller, String msg) {
        error(caller, msg, System.err);
    }

    public static void error(Object caller, String msg, PrintStream ps) {
        ps.println( "[ERROR]\t" + caller.getClass().getSimpleName() + "\t" + timestamp() + " - " + msg);
        ps.flush();
    }

    public static void debug(Object caller, String msg) {
        if (verbose) {
            debug(caller, msg, System.err);
        }
    }

    public static void debug(Object caller, String msg, PrintStream ps) {
        if (verbose) {
            ps.println("[DEBUG]\t" + caller.getClass().getSimpleName() + "\t" + timestamp() + " - " + msg);
            ps.flush();
        }
    }

    public static String timestamp() {
        return timestamp(DATE_FORMAT);
    }

    public static String timestamp(SimpleDateFormat dataFormat) {
        return dataFormat.format(Calendar.getInstance().getTime());
    }
}
