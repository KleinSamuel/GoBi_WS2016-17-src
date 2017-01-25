package bamfiles;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import assignment_3.Task_3_Gene_Counts;
import debugStuff.DebugMessageFactory;
import genomeAnnotation.Gene;
import genomeAnnotation.GenomeAnnotation;
import genomeAnnotation.Transcript;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import io.ConfigReader;
import javafx.util.Pair;
import reader.GTFParser;

public class BAMFileReader {

	// readId, waitingRead
	private HashMap<String, SAMRecord> waitingRecords;
	private String bamFile;
	public static GenomeAnnotation ga;
	
	private Counter counter;
	private HashMap<String, Task_3_Gene_Counts> map;

	public BAMFileReader(String bamPath, Counter counter, HashMap<String, Task_3_Gene_Counts> map) {
		bamFile = bamPath;
		waitingRecords = new HashMap<>();
		this.counter = counter;
		this.map = map;
		ga = GTFParser.readGtfFile("h.ens.75", "/home/proj/biosoft/praktikum/genprakt-ws16/gtf/Homo_sapiens.GRCh37.75.gtf");
		
	}

	public void readBAMFile() {
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader sr = SamReaderFactory.makeDefault().open(new File(bamFile));
		Iterator<SAMRecord> it = sr.iterator();
		SAMRecord sam = null, possibleMate = null;
		// reads are sorted by start --> so if new chromosome clear map
		String chromId = null;
		ReadPair rp = null;
		int validRecords = 0, validPairs = 0, invalidRecords = 0, nonValidPairs = 0, checkedRecords = 0;
		
		while (it.hasNext()) {
			
			sam = it.next();
			checkedRecords++;
			
			if (checkedRecords % 100000 == 0) {
				DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE,
						"checkedRecords: " + checkedRecords + "\tvalidRecords: " + validRecords + "\tinvalidRecords: "
								+ invalidRecords + "\tvalidPairs: " + validPairs + "\tnonValidPairs: " + nonValidPairs);
			}
			
			if (validRecord(sam)) {
				
				validRecords++;
				
				// check if new chromosome
				if (chromId == null) {
					chromId = sam.getReferenceName();
				} else {
					if (!sam.getReferenceName().equals(chromId)) {
						waitingRecords = new HashMap<>();
					}
				}
				
				// look for waiting record in map
				possibleMate = waitingRecords.get(sam.getReadName());
				
				if (possibleMate == null) {
					waitingRecords.put(sam.getReadName(), sam);
				} else {
					// check if valid pair --> but what to do if not --> both
					// reads valid?? possible??
					rp = validPair(sam, possibleMate);
					if (rp == null) {
						rp = validPair(possibleMate, sam);
					}
					if (rp == null) {
						nonValidPairs++;
						continue;
					} else {
						waitingRecords.remove(sam.getReadName());
						validPairs++;
						
						/* task 2 */
						if(counter != null){
							
							this.counter.addNRP();
							
							if(rp.getMatchedGenes().size() == 1){
								this.counter.addMappedNRP();
							}
							if(rp.getMatchedGenes().size() > 1){
								this.counter.addMultimappedNRP();
							}
							
							if(rp.isTranscriptomic()){
								this.counter.addTranscriptomicNRP();
							}
							if(rp.isMerged()){
								this.counter.addMergedTranscriptNRP();
							}
							if(rp.isIntronic()){
								this.counter.addIntronicNRP();
							}
							if(rp.isAntisense()){
								this.counter.addAntisenseNRP();
							}
							if(rp.isIntergenic()){
								this.counter.addIntergenicNRP();
							}
							
						}
						
					}
				}
			} else {
				invalidRecords++;
			}
		}
		
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE,
				"checkedRecords: " + checkedRecords + "\tvalidRecords: " + validRecords + "\tinvalidRecords: "
						+ invalidRecords + "\tvalidPairs: " + validPairs + "\tnonValidPairs: " + nonValidPairs);
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "Finished reading");

	}
	
	public Counter getCounter(){
		return this.counter;
	}

	public boolean validRecord(SAMRecord sam) {
		return (!sam.getReadUnmappedFlag() && !sam.getMateUnmappedFlag() && !sam.getNotPrimaryAlignmentFlag()
				&& sam.getReferenceName().equals(sam.getMateReferenceName())
				&& sam.getReadNegativeStrandFlag() != sam.getMateNegativeStrandFlag());
	}

	public ReadPair validPair(SAMRecord first, SAMRecord second) {
		if (!first.getReferenceName().equals(second.getReferenceName()))
			return null;
		if (first.getFirstOfPairFlag() && second.getSecondOfPairFlag()
				&& first.getAlignmentStart() == second.getMateAlignmentStart()
				&& first.getMateAlignmentStart() == second.getAlignmentStart()) {
			return new ReadPair(first, second, false);
		}
		return null;
	}

}