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
package com.gkano.bioinfo.var;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.RandomAccessFile;
import java.net.URISyntaxException;
import java.net.URL;
import java.security.SecureRandom;
import java.text.DateFormat;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.PrimitiveIterator.OfInt;
import java.util.zip.GZIPInputStream;

public class GeneralTools {

    public static final DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
    public static final DecimalFormat decimalFormat = new DecimalFormat("#.######", new DecimalFormatSymbols(Locale.US));

    private static OfInt prng;

    static {
        try {
            //prng = SecureRandom.getInstance("SHA1PRNG");
            SecureRandom random = new SecureRandom();
            prng = random.ints().distinct().iterator();
        } catch (Exception e) {
            Logger.error(GeneralTools.class, e.getMessage());
        }
    }

    private static GeneralTools instance;

    private GeneralTools() {
    }

    public static GeneralTools getInstance() {
        if (instance == null) {
            instance = new GeneralTools();
        }
        return instance;
    }

    public static int getRandomInteger() {
        return prng.nextInt();
    }

    public static PrintStream getPrintStreamOrExit(String outputFile, Object context) {
        try {
            return (outputFile == null) ? System.out : new PrintStream(outputFile);
        } catch (FileNotFoundException | SecurityException e) {
            Logger.error(context, "Cannot write to " + outputFile);
            System.exit(1);
            return null;
        }
    }

    public static String time() {
        Calendar cal = Calendar.getInstance();
        return dateFormat.format(cal.getTime());
    }

    public static String RAMInfo(Runtime runtime) {
        double gb = (double) 1024 * 1024 * 1024;
        NumberFormat format = NumberFormat.getInstance();
        StringBuilder sb = new StringBuilder();

        sb.append("\n##### Heap utilization statistics [GB] #####\n");
        sb.append("Max Memory=").append(format.format(runtime.maxMemory() / gb)).append("\n");
        sb.append("Current Total Memory=").append(format.format(runtime.totalMemory() / gb)).append("\n");
        sb.append("Current Used Memory=").append(format.format((runtime.totalMemory() - runtime.freeMemory()) / gb)).append("\n");
        sb.append("Current Free Memory=").append(format.format(runtime.freeMemory() / gb)).append("\n");
        sb.append("############################################\n");

        return sb.toString();
    }

    final static public byte[][] concatSeqQual(final ArrayList<byte[]> arrays) {
        byte[] all = concat(arrays);
        int size = all.length;
        int sizeSeq = 0;
        int sizeQual = 0;
        for (byte b : all) {
            if ('+' == (char) (b & 0xff)) {
                sizeQual = size - sizeSeq - 1;
                break;
            } else {
                sizeSeq++;
            }
        }

        byte[][] ret = new byte[2][];
        if (sizeSeq > 0) {
            byte[] destSeq = new byte[sizeSeq];
            System.arraycopy(all, 0, destSeq, 0, sizeSeq);
            ret[0] = destSeq;
        }
        if (sizeQual > 0) {
            byte[] destQual = new byte[sizeQual];
            System.arraycopy(all, sizeSeq + 1, destQual, 0, sizeQual);
            ret[1] = destQual;
        }

        return ret;
    }

    final static public byte[] concat(final ArrayList<byte[]> arrays) {
        int size = 0;
        for (byte[] a : arrays) {
            size += a.length;
        }

        byte[] dest = new byte[size];

        int destPos = 0;
        for (int i = 0; i < arrays.size(); i++) {
            if (i > 0) {
                destPos += arrays.get(i - 1).length;
            }
            int length = arrays.get(i).length;
            System.arraycopy(arrays.get(i), 0, dest, destPos, length);
        }

        return dest;
    }

    public static int round(double d) {
        double dAbs = Math.abs(d);
        int i = (int) dAbs;
        double result = dAbs - (double) i;
        if (result < 0.5) {
            return d < 0 ? -i : i;
        } else {
            return d < 0 ? -(i + 1) : i + 1;
        }
    }

    public static char[] bytesToStringUTFCustom(byte[] bytes) {
        char[] buffer = new char[bytes.length >> 1];
        for (int i = 0; i < buffer.length; i++) {
            int bpos = i << 1;
            char c = (char) (((bytes[bpos] & 0x00FF) << 8) + (bytes[bpos + 1] & 0x00FF));
            buffer[i] = c;
        }
        return buffer;
        //return new String(buffer);
    }

    public static byte[] stringToBytesUTFCustom(char[] str) {
        byte[] b = new byte[str.length << 1];
        for (int i = 0; i < str.length; i++) {
            char strChar = str[i];
            int bpos = i << 1;
            b[bpos] = (byte) ((strChar & 0xFF00) >> 8);
            b[bpos + 1] = (byte) (strChar & 0x00FF);
        }
        return b;
    }

    public static List<String> directoriesToFiles(List<String> inputFastaFileNames) {
        ArrayList<String> outputFastaFileNames = new ArrayList<>();
        try {
            for (String s : inputFastaFileNames) {
                File f = new File(s);
                if (!f.isDirectory()) {
                    outputFastaFileNames.add(f.getCanonicalPath());
                } else {
                    for (File ff : f.listFiles()) {
                        outputFastaFileNames.add(ff.getCanonicalPath());
                    }
                }
            }
        } catch (IOException e) {
            Logger.error(GeneralTools.class, e.getMessage());
        }
        return outputFastaFileNames;
    }

    public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
        List<Map.Entry<K, V>> list = new LinkedList<>(map.entrySet());
        Collections.sort(list, (Map.Entry<K, V> o1, Map.Entry<K, V> o2) -> {
            int compare = (o1.getValue()).compareTo(o2.getValue());
            if (compare == 0) {
                Long key1 = (Long) o1.getKey();
                Long key2 = (Long) o2.getKey();
                return +1 * (key1).compareTo(key2);
            } else {
                return -1 * compare;
            }
        });

        Map<K, V> result = new LinkedHashMap<>();
        for (Map.Entry<K, V> entry : list) {
            result.put(entry.getKey(), entry.getValue());
        }
        return result;
    }

    public static ArrayList<String> generateKmerVocabulary(int k) {
        ArrayList<String> vocabulary = new ArrayList<>();

        char[] chars = "ATCG".toCharArray();
        int len = k;
        iterate(chars, len, new char[len], 0, vocabulary);

        return vocabulary;
    }

    private static void iterate(char[] chars, int len, char[] build, int pos, ArrayList<String> vocabulary) {
        if (pos == len) {
            String word = new String(build);
            vocabulary.add(word);
            return;
        }

        for (int i = 0; i < chars.length; i++) {
            build[pos] = chars[i];
            iterate(chars, len, build, pos + 1, vocabulary);
        }
    }

    /*
	public static String getRandomSequence(int length){
		return RandomStringUtils.random(length, "ATCG");
	}
     */
    public static final String reverseComplement(String sequence) {
        StringBuilder sb = new StringBuilder();
        for (int i = sequence.length() - 1; i >= 0; i--) {
            switch (sequence.charAt(i)) {
                case 'A':
                case 'a':
                    sb.append("T");
                    break;

                case 'T':
                case 't':
                    sb.append("A");
                    break;

                case 'C':
                case 'c':
                    sb.append("G");
                    break;

                case 'G':
                case 'g':
                    sb.append("C");
                    break;

                default:
					;
                    break;
            }
        }
        return sb.toString();
    }

    @SuppressWarnings("rawtypes")
    public static long classBuildTimeMillis(Class c) {
        URL resource = c.getResource(c.getSimpleName() + ".class");
        if (resource == null) {
            System.out.println("Failed to find class file for class: " + c.getName());
            return 0L;
        }

        switch (resource.getProtocol()) {
            case "file":
                try {
                    return new File(resource.toURI()).lastModified();
                } catch (URISyntaxException e) {
                    Logger.error(GeneralTools.class, e.getMessage());
                    return 0L;
                }
            case "jar":
                String path = resource.getPath();
                return new File(path.substring(5, path.indexOf("!"))).lastModified();
            default:
                System.out.println("Unhandled url protocol: " + resource.getProtocol() + " for class: " + c.getName() + " resource: " + resource.toString());
                return 0L;
        }
    }

    public static <T extends Comparable<? super T>> List<T> asSortedList(Collection<T> c) {
        List<T> list = new ArrayList<>(c);
        Collections.sort(list);
        return list;
    }

    public static String[] reorderLabels(String[] original, String treeString) {
        String[] labelsReordered = new String[original.length];

        String[] data = treeString.split(",");
        for (int i = 0; i < data.length; i++) {
            String label = data[i].trim();

            if (label.contains("(")) {
                label = label.substring(label.lastIndexOf("(") + 1);
            }
            if (label.contains(")")) {
                label = label.substring(0, label.indexOf(")"));
            }
            if (label.contains(":")) {
                label = label.substring(0, label.indexOf(":"));
            }

            //System.out.println("i="+i+"\tlabel="+label);
            labelsReordered[i] = label;
        }

        return labelsReordered;
    }

    public static double[][] reorderDistances(double[][] distancesOriginal, String[] labelsOriginal, String[] labelsReordered) {
        int size = distancesOriginal.length;
        double[][] distancesReordered = new double[size][size];

        int[] indexToReorder = new int[size];
        for (int i = 0; i < size; i++) {
            String labelOriginal = labelsOriginal[i].trim();
            for (int j = 0; j < size; j++) {
                String labelReordered = labelsReordered[j].trim();
                if (labelOriginal.equalsIgnoreCase(labelReordered)) {
                    indexToReorder[i] = j;
                    break;
                }
            }
        }

        //reorder
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                //distancesReordered[i][j] = distancesOriginal[indexToReorder[i]][indexToReorder[j]];
                distancesReordered[indexToReorder[i]][indexToReorder[j]] = distancesOriginal[i][j];
            }
        }

        return distancesReordered;
    }

    public static Object[] readDistancesSamples(String input) {
        try {
            Object[] data;
            try (FileInputStream fis = new FileInputStream(input)) {
                data = readDistancesSamples(fis);
            }
            return data;
        } catch (IOException e) {
            Logger.error(GeneralTools.class, e.getMessage());
            return null;
        }
    }

    public static Object[] readDistancesSamples(InputStream input) {
        try {
            double[][] distances;
            String[] sampleNames;
            try (BufferedReader br = new BufferedReader(new InputStreamReader(input))) {
                String line = br.readLine();
                System.err.println(line);
                int numOfSamples = Integer.parseInt(line.split("[\\s,\\t]+")[0]);
                distances = new double[numOfSamples][numOfSamples];
                sampleNames = new String[numOfSamples];
                int i = 0;
                while ((line = br.readLine()) != null) {
                    String[] data = line.split("[\\s,\\t]+");
                    sampleNames[i] = data[0].trim();
                    for (int j = 0; j < numOfSamples; j++) {
                        try {
                            distances[i][j] = Double.parseDouble(data[j + 1]);
                        } catch (NumberFormatException e) {
                            distances[i][j] = 1.0;
                        }
                    }
                    i++;
                }
            }

            Object[] ret = new Object[2];
            ret[0] = distances;
            ret[1] = sampleNames;

            return ret;
        } catch (IOException | NumberFormatException e) {
            Logger.error(GeneralTools.class, e.getMessage());
            return null;
        }
    }

    public static boolean isGZipped(File f) {
        int magic = 0;
        try {
            try (RandomAccessFile raf = new RandomAccessFile(f, "r")) {
                magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
            }
        } catch (IOException e) {
            Logger.error(GeneralTools.class, e.getMessage());
        }
        return magic == GZIPInputStream.GZIP_MAGIC;
    }

    public static int getAdaptiveVariantStep(int count) {
        if (count < 10_000) {
            return 100;
        } else if (count < 100_000) {
            return 1_000;
        } else if (count < 1_000_000) {
            return 10_000;
        } else {
            return 100_000;
        }
    }

    public static int getAdaptiveSampleStep(int count) {
        if (count < 100) {
            return 10;
        } else if (count < 1_000) {
            return 100;
        } else {
            return 1_000;
        }
    }
}
