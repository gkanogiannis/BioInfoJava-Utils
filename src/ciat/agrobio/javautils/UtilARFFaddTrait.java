package ciat.agrobio.javautils;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;

import ciat.agrobio.core.GeneralTools;

@Parameters(commandDescription = "ARFFaddTrait")
public class UtilARFFaddTrait {
	
	private static UtilARFFaddTrait instance = new UtilARFFaddTrait();
	
	private UtilARFFaddTrait() {}
	
	public static UtilARFFaddTrait getInstance() {return instance;}
	
	public static String getUtilName() {
		return "ARFFaddTrait";
	}
	
	@Parameter(names = "--help", help = true)
	private boolean help;
	
	@Parameter(names={"--inputARFF","-i"}, required=true)
	private String inputFileName;
	
	@Parameter(names={"--inputTrait","-t"}, required=true)
	private String traitFileName;
	
	@Parameter(names={"--traitAttributeString","-n"}, required=true, description="classExample numeric")
	private String traitString;
	
	public void go() {
		try {
			InputStream fis1 = Files.newInputStream(Paths.get(inputFileName));
			InputStream fis2 = Files.newInputStream(Paths.get(traitFileName));
			BufferedReader br1,br2;
			if(inputFileName.endsWith(".gz")) {
				br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis1),"UTF-8"));
			}
			else {
				br1 = new BufferedReader(new InputStreamReader(fis1, "UTF-8"));
			}
			if(traitFileName.endsWith(".gz")) {
				br2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(fis2),"UTF-8"));
			}
			else {
				br2 = new BufferedReader(new InputStreamReader(fis2, "UTF-8"));
			}
		
			String line1; 
			int sampleCounter = 0;
			while((line1=br1.readLine()) != null) {
				if(line1.startsWith("@data")) {
					System.out.println("@attribute\t" + traitString + "\n");
					System.out.println("@data");
					break;
				}
				else {
					System.out.println(line1);
				}
			}
			String line2; 
			while((line1=br1.readLine()) != null && (line2=br2.readLine()) != null) {
				try {
					if(line2==null || line2.isEmpty()) System.out.println(line1+",?");
					else {
						System.out.println(line1+","+line2);
					}
				} 
				catch (Exception e) {
					System.out.println(line1+",?");
				}
				if(++sampleCounter % 10 == 0) System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			}
			
			System.err.println(GeneralTools.time()+" Samples Processed : \t"+sampleCounter);
			br1.close();
			br2.close();
		} 
		catch (Exception e) {
			e.printStackTrace();
		}
	}
}

