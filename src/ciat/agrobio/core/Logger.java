package ciat.agrobio.core;

import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class Logger {

    private static final SimpleDateFormat DATE_FORMAT = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private static boolean verbose = false;

    public static void setVerbose(boolean v) {
        verbose = v;
    }

    public static void info(String msg) {
        if (verbose) info(msg, System.err);
    }
    public static void info(String msg, PrintStream ps) {
        if (verbose) ps.println("[INFO] " + timestamp() + " - " + msg);
    }
    
    public static void infoCarret(String msg) {
        if (verbose) info(msg, System.err);
    }
    public static void infoCarret(String msg, PrintStream ps) {
        if (verbose) {
            ps.print("\r[INFO] " + timestamp() + " - " + msg);
            ps.flush();
        }
    }

    public static void warn(String msg) {
        warn(msg, System.err);
    }
    public static void warn(String msg, PrintStream ps) {
        ps.println("[WARN] " + timestamp() + " - " + msg);
    }

    public static void error(String msg) {
        error(msg, System.err);
    }
    public static void error(String msg, PrintStream ps) {
        ps.println("[ERROR] " + timestamp() + " - " + msg);
    }

    public static void debug(String msg) {
        if (verbose) debug(msg, System.err);
    }
    public static void debug(String msg, PrintStream ps) {
        if (verbose) ps.println("[DEBUG] " + timestamp() + " - " + msg);
    }

    public static String timestamp() {
        return timestamp(DATE_FORMAT);
    }
    public static String timestamp(SimpleDateFormat dataFormat) {
        return dataFormat.format(Calendar.getInstance().getTime());
    }
}
