package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Iterator;

import assignment_3.Task_2_Counts;
import assignment_3.Task_3_Gene_Counts_Container;
import assignment_3.Task_4_ExtendedFeatureStatistic;
import assignment_3.splicing_evidence.SplicingEvidences;
import bamfiles.PCRIndexStructure;
import bamfiles.ReadPair;
import debugStuff.DebugMessageFactory;
import genomeAnnotation.GenomeAnnotation;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import io.ConfigReader;
import io.HeadedFileReader.LineObject;
import reader.GTFParser;

public class Anders {

	// readId, waitingRead
	private HashMap<String, SAMRecord> waitingRecords;
	private String bamFile;
	public static GenomeAnnotation ga;
	private String ebv_mapper;
	private LineObject mapping;

	public Anders(String bamPath, String gtfPath, LineObject mapping) {
		bamFile = bamPath;
		waitingRecords = new HashMap<>();
		ga = GTFParser.readGtfFile("h.ens.75", gtfPath);
		ebv_mapper = mapping.getValue("condition").replace(":", "/") + "/rep" + mapping.getValue("replicate") + "/"
				+ mapping.getValue("mapper") + "/" + mapping.getValue("mapper");
		File f = new File(ConfigReader.readConfig().get("output_directory") + "/"
				+ ebv_mapper.substring(0, ebv_mapper.lastIndexOf("/") + 1));
		f.mkdirs();
		f = null;
		this.mapping = mapping;
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
		int splitInconsistent = 0, wrong = 0;
		PCRIndexStructure pcrIndices = new PCRIndexStructure();

		SplicingEvidences se = new SplicingEvidences(mapping.getValue("mapper"), ebv_mapper);

		Task_2_Counts featureCounts = new Task_2_Counts(Integer.parseInt(mapping.getValue("nreads")),
				mapping.getValue("mapper"), ebv_mapper, new Task_4_ExtendedFeatureStatistic());

		Task_3_Gene_Counts_Container geneCounts = new Task_3_Gene_Counts_Container(ebv_mapper, featureCounts);

		BufferedWriter bw = null;
		try {
			bw = new BufferedWriter(new FileWriter(
					new File(ConfigReader.readConfig().get("output_directory") + ebv_mapper + ".annot")));
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		while (it.hasNext()) {
			sam = it.next();
			checkedRecords++;
			if (mappedRecordAndMate(sam)) {
				featureCounts.increaseMapped();
			}
			if (multimappedRecordAndMate(sam)) {
				featureCounts.increaseMultimapped();
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
						if (rp.getSplitCount() < 0) {
							splitInconsistent++;
							try {
								bw.write(rp.getReadName() + "\tsplit-inconsistent:true\n");
							} catch (Exception e) {
								e.printStackTrace();
								System.exit(1);
							}
						} else {
							int pcrIndex = pcrIndices.addReadPair(rp, bw);
							rp.setPCRindex(pcrIndex);
							try {
								bw.write(rp.getReadName() + "\t" + rp.getAttributeXX() + "\n");

								// if (!compareXX(rp)) {
								// bw.write(rp.getAttributeXX() + "\n");
								// bw.write(rp.XX + "\n");
								// for (Gene g :
								// ga.getChromosome(rp.getReferenceName())
								// .getSpecificStrandGenes(rp.isOnNegativeStrand())
								// .getIntervalsIntersecting(rp.getGenomicRegion().getStart(),
								// rp.getGenomicRegion().getStop(), new
								// LinkedList<>())) {
								// bw.write(rp.getReadName() + "\t" +
								// rp.getGenomicRegion().toString() + "\t"
								// + rp.isOnNegativeStrand() + "\n");
								// bw.write(g.getId() + ": " + g.getStart() +
								// "-" + g.getStop() + "\t"
								// + g.isOnNegativeStrand() + "|");
								// bw.write("\n");
								// }
								// }
							} catch (Exception e) {
								e.printStackTrace();
								System.exit(1);
							}
							geneCounts.addReadPair(rp);
							se.addIntrons(rp);
						}
						featureCounts.addReadPair(rp);
					}
				}
			} else {
				invalidRecords++;
			}
			if (checkedRecords % 100000 == 0) {
				DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE,
						"checkedRecords: " + checkedRecords + "\tvalidRecords: " + validRecords + "\tinvalidRecords: "
								+ invalidRecords + "\tvalidPairs: " + validPairs + "\tnonValidPairs: " + nonValidPairs
								+ "\tinconsistent: " + splitInconsistent + "\twrong: " + wrong);
			}
		}
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE,
				"checkedRecords: " + checkedRecords + "\tvalidRecords: " + validRecords + "\tinvalidRecords: "
						+ invalidRecords + "\tvalidPairs: " + validPairs + "\tnonValidPairs: " + nonValidPairs
						+ "\tinconsistent: " + splitInconsistent + "\twrong: " + wrong);
		DebugMessageFactory.printInfoDebugMessage(ConfigReader.DEBUG_MODE, "Finished reading");
		try {
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		se.close();
		se.plot();
		featureCounts.barPlot();
		geneCounts.plot();
		geneCounts.writeGeneCountsTSV();
	}

	public boolean compareXX(ReadPair rp) {
		return rp.compareXX();
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

	public boolean mappedRecordAndMate(SAMRecord sam) {
		if (!sam.getReadUnmappedFlag() && !sam.getMateUnmappedFlag()) {
			return true;
		}
		return false;
	}

	public boolean multimappedRecordAndMate(SAMRecord sam) {
		if (!sam.getReadUnmappedFlag() && !sam.getMateUnmappedFlag()) {
			if (sam.isSecondaryOrSupplementary()) {
				return true;
			}
		}
		return false;
	}

}
