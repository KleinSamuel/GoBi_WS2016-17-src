package bamfiles;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.Vector;

import assignment_3.BAMFilePaths;
import assignment_3.Gene_Counts;
import io.ConfigHelper;
import io.ExternalFileWriter;
import javafx.util.Pair;
import plotting.BarPlot;

public class BAMStatistics {

	public HashMap<String, Gene_Counts> countMap;
	public static HashMap<String, HashMap<String, Integer>> transCountMap;
	public static HashMap<String, HashMap<String, Integer>> intronCountMap;
	public static HashMap<String, Integer> mergedTrLength;
	
	public void task2(){
		
		ArrayList<String> paths = BAMFilePaths.getPathList();
		
		for(String s : paths){
			
			s = "/home/proj/biosoft/praktikum/genprakt-ws16/assignment/a3/data/debug_bams/contextmap.bam";
			
			System.out.println(s);
			BAMFileReader br = new BAMFileReader(s, new Counter(), null);
			br.readBAMFile();
			Counter ct = br.getCounter();
			System.out.println(ct.toString());
			
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
			
//			String name = s.substring(s.indexOf("EBNA2/"));
			
			BarPlot bp = new BarPlot(pair, "Test", "", "Amount", false);
			bp.filename = "TestTask2";
			bp.plot();
			
			try {
				Thread.sleep(1000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			System.exit(0);
		}
		
	}
	
	public static int allNRP;
	public static int mergedTrNRP;
	
	public void task3(){
		
		ArrayList<String> paths = BAMFilePaths.getPathList();
		
		
		for(String s : paths){
			
			this.countMap = new HashMap<>();
			this.transCountMap = new HashMap<>();
			
			TreeMap<Double, Integer> geneTrans = new TreeMap<>();
			TreeMap<Double, Integer> geneIntron = new TreeMap<>();
			
			double gene_tr_FPKM = 0.0;
			double gene_intron_FPKM = 0.0;
			
			s = "/home/proj/biosoft/praktikum/genprakt-ws16/assignment/a3/data/debug_bams/contextmap.bam";
			
			System.out.println(s);
			BAMFileReader br = new BAMFileReader(s, null, this.countMap);
			br.readBAMFile();
			
//			String name = s.substring(s.indexOf("EBNA2/"));
			
			ExternalFileWriter ew = new ExternalFileWriter();
			ew.openWriter(new ConfigHelper().getDefaultOutputPath()+"counts.tsv");
			
			for(Entry<String, Gene_Counts> entry : countMap.entrySet()){
				
				ew.writeToWriter(entry.getKey()+"\t");
				ew.writeToWriter(entry.getValue().getNRPwithinGenReg()+"\t");
				ew.writeToWriter(entry.getValue().getNRPwithinGenRegWitoutOtherGenReg()+"\t");
				ew.writeToWriter(entry.getValue().getNRPintronic()+"\t");
				ew.writeToWriter(entry.getValue().getNRPtranscriptomic()+"\t");
				ew.writeToWriter(entry.getValue().getNRPmergedTr()+"\t");
				
				HashMap<String, Integer> tmp = transCountMap.get(entry.getKey());
				Integer max = 0;
				
				if(tmp != null){
					for(Integer i : tmp.values()){
						max = Math.max(max, i);
					}
				}
				
				ew.writeToWriter(max+"\n");
				
				/* gene-tr-FPKM plot */
//				gene_tr_FPKM = (Math.pow(10, 9) * entry.getValue().getNRPtranscriptomic())/(allNRP*mergedTrLength.get(entry.getKey()));
//				
//				if(geneTrans.containsKey(gene_tr_FPKM)){
//					geneTrans.put(gene_tr_FPKM, geneTrans.get(gene_tr_FPKM)+1);
//				}else{
//					geneTrans.put(gene_tr_FPKM, 1);
//				}
				
				/* gene-intron-FPKM plot */
				if(mergedTrLength.containsKey(entry.getKey())){
					gene_intron_FPKM = (Math.pow(10, 9)*entry.getValue().getNRPintronic())/(allNRP*mergedTrLength.get(entry.getKey()));
					
					if(geneIntron.containsKey(gene_intron_FPKM)){
						geneIntron.put(gene_intron_FPKM, geneIntron.get(gene_intron_FPKM)+1);
					}else{
						geneIntron.put(gene_intron_FPKM, 1);
					}
				}
				
			}
			
			for(Entry<Double, Integer> entry : geneIntron.entrySet()){
				System.out.println(entry.getKey()+"\t:\t"+entry.getValue());
			}
			
			ew.closeWriter();
			System.exit(0);
		}
		
	}
	
	public static void main(String[] args) {
		
		BAMStatistics bs = new BAMStatistics();
		bs.task3();
		
	}
	
	
}
