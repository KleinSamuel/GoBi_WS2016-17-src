package bamfiles;

import java.util.ArrayList;
import java.util.Vector;

import assignment_3.BAMFilePaths;
import javafx.util.Pair;
import plotting.BarPlot;

public class BAMStatistics {

	
	public void readBamFile(){
		
		ArrayList<String> paths = BAMFilePaths.getPathList();
		
		for(String s : paths){
			System.out.println(s);
			BAMFileReader br = new BAMFileReader(s, new Counter());
			br.readBAMFile();
			Counter ct = br.getCounter();
			
			Vector<Object> key = new Vector<>();
			Vector<Object> value = new Vector<>();
			
			key.add(ct.getNRP());
			value.add("Number Read Pairs");
			key.add(ct.getMappedNRP());
			value.add("Mapped NRP");
			key.add(ct.getMultimappedNRP());
			value.add("Multimapped NRP");
			key.add(ct.getTranscriptomicNRP());
			value.add("Transcriptomic NRP");
			key.add(ct.getMergedTranscriptNRP());
			value.add("Merged_tr NRP");
			key.add(ct.getIntronicNRP());
			value.add("Intronic NRP");
			key.add(ct.getAntisenseNRP());
			value.add("Antisense NRP");
			key.add(ct.getIntergenicNRP());
			value.add("Intergenic NRP");
			
			Pair<Vector<Object>, Vector<Object>> pair = new Pair<Vector<Object>, Vector<Object>>(key, value);
			
			String name = s.substring(s.indexOf("EBNA2/"));
			
			BarPlot bp = new BarPlot(pair, name, "", "Amount", false);
			bp.filename = name;
			bp.plot();
		}
		
	}
	
	public static void main(String[] args) {
		
		BAMStatistics bs = new BAMStatistics();
		bs.readBamFile();
		
	}
	
	
}
