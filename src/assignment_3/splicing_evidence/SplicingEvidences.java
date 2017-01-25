package assignment_3.splicing_evidence;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.TreeMap;

import bamfiles.ReadPair;
import io.ConfigReader;
import util.Interval;

public class SplicingEvidences {

	// nur bei merged und intronic checken

	private TreeMap<Interval, SplicingEvidence> posStrandIntrons, negStrandIntrons;
	private String currentReferenceName = "NoRefThereYet";
	private BufferedWriter bw;
	private Splicing_Counts counts;
	private String ebv_mapper, mapper;

	/**
	 * has to be closed when adding all readpairs is finished
	 * 
	 * @param ga
	 * @param mappername
	 */
	public SplicingEvidences(String mapper, String ebv_mappername) {
		posStrandIntrons = new TreeMap<>();
		negStrandIntrons = new TreeMap<>();
		counts = new Splicing_Counts();
		this.mapper = mapper;
		this.ebv_mapper = ebv_mappername;
		try {
			bw = new BufferedWriter(new FileWriter(new File(
					ConfigReader.readConfig().get("output_directory") + "/" + ebv_mappername + "_splitinfos.tsv")));
			writeHeader();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void addIntrons(ReadPair rp) {
		if (rp.getSplitCount() < 1) {
			return;
		}
		if (!rp.getReferenceName().equals(currentReferenceName)) {
			for (SplicingEvidence s : negStrandIntrons.values()) {
				writeSplicingEvidenceToFile(s);
			}
			for (SplicingEvidence s : posStrandIntrons.values()) {
				writeSplicingEvidenceToFile(s);
			}
			posStrandIntrons = new TreeMap<>();
			negStrandIntrons = new TreeMap<>();
			currentReferenceName = rp.getReferenceName();
		}
		if (negStrandIntrons.size() > 200000
				&& rp.getGenomicRegion().getStart() > negStrandIntrons.lastKey().getStart() / 2) {
			Iterator<Interval> iter = negStrandIntrons.keySet().iterator();
			SplicingEvidence se = null;

			while (iter.hasNext()) {
				se = negStrandIntrons.get(iter.next());
				if (se.getStart() < rp.getGenomicRegion().getStart() / 2) {
					writeSplicingEvidenceToFile(se);
					iter.remove();
				}
			}
		}
		if (posStrandIntrons.size() > 200000
				&& rp.getGenomicRegion().getStart() > posStrandIntrons.lastKey().getStart() / 2) {
			Iterator<Interval> iter = posStrandIntrons.keySet().iterator();
			SplicingEvidence se = null;

			while (iter.hasNext()) {
				se = posStrandIntrons.get(iter.next());
				if (se.getStart() < rp.getGenomicRegion().getStart() / 2) {
					writeSplicingEvidenceToFile(se);
					iter.remove();
				}
			}
		}
		LinkedList<Interval> intronsFor = rp.getForwardIntrons(), intronsRev = rp.getReverseIntrons();
		for (Interval i : intronsFor) {
			addRP(rp, i);
		}
		for (Interval i : intronsRev) {
			addRP(rp, i);
		}
	}

	public void addRP(ReadPair rp, Interval i) {
		SplicingEvidence next = null;
		if (rp.isOnNegativeStrand()) {
			next = negStrandIntrons.get(i);
		} else {
			next = posStrandIntrons.get(i);
		}
		if (next == null) {
			next = new SplicingEvidence(i);
			if (rp.isOnNegativeStrand()) {
				negStrandIntrons.put(i, next);
			} else {
				posStrandIntrons.put(i, next);
			}
		}
		next.addReadPair(rp);
	}

	public void writeHeader() {
		try {
			bw.write(
					"chr\tstart\tend\tstrand\tgene_id\ttranscript_ids\tannot\tcount\tucounts\tcounts_0\tucounts_0\tcounts_1\tucounts_1\tcounts_2\tucounts_2\n");
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void writeSplicingEvidenceToFile(SplicingEvidence se) {
		try {
			bw.write(se.toString() + "\n");
			counts.addSplicingEvidence(se);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void close() {
		for (SplicingEvidence s : negStrandIntrons.values()) {
			writeSplicingEvidenceToFile(s);
		}
		for (SplicingEvidence s : posStrandIntrons.values()) {
			writeSplicingEvidenceToFile(s);
		}
		negStrandIntrons = null;
		posStrandIntrons = null;
		try {
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void plot() {
		counts.plot(mapper, ebv_mapper);
	}

}
