package assignment_4;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Map.Entry;
import java.util.TreeMap;

public class EBDataCollector {

	private TreeMap<String, int[]> geneCounts;
	private TreeMap<String, TreeMap<String, int[]>> transcriptCounts;

	public static void main(String[] args) {
		new EBDataCollector();
	}

	public EBDataCollector() {
		geneCounts = new TreeMap<>();
		transcriptCounts = new TreeMap<>();
		ArrayList<String> files = new ArrayList<>(6);
		files.add("D:/Dennis/Uni/GoBi/A4/diff_simulation/f1/r1/read.mappinginfo");
		files.add("D:/Dennis/Uni/GoBi/A4/diff_simulation/f1/r2/read.mappinginfo");
		files.add("D:/Dennis/Uni/GoBi/A4/diff_simulation/f1/r3/read.mappinginfo");
		files.add("D:/Dennis/Uni/GoBi/A4/diff_simulation/f2/r1/read.mappinginfo");
		files.add("D:/Dennis/Uni/GoBi/A4/diff_simulation/f2/r2/read.mappinginfo");
		files.add("D:/Dennis/Uni/GoBi/A4/diff_simulation/f2/r3/read.mappinginfo");
		for (int i = 0; i < files.size(); i++) {
			readCountFile(i, files.get(i));
		}
		writeFeatures();
		writeExpression();
		writeTranscriptFeatures();
		writeTranscriptExpression();

	}

	public void readCountFile(int index, String file) {
		try {
			String line = null;
			String[] split = null;
			BufferedReader br = new BufferedReader(new FileReader(new File(file)));
			br.readLine();
			while ((line = br.readLine()) != null) {
				split = line.split("\t");
				int[] rcounts = geneCounts.get(split[0]);
				if (rcounts == null) {
					rcounts = new int[] { 0, 0, 0, 0, 0, 0 };
					geneCounts.put(split[0], rcounts);
				}
				rcounts[index] += Integer.parseInt(split[2]);

				// trCounts
				TreeMap<String, int[]> trCounts = transcriptCounts.get(split[0]);
				int[] trCountArr = null;
				if (trCounts == null) {
					trCounts = new TreeMap<>();
					trCountArr = new int[] { 0, 0, 0, 0, 0, 0 };
					trCounts.put(split[1], trCountArr);
					transcriptCounts.put(split[0], trCounts);
				} else {
					trCountArr = trCounts.get(split[1]);
					if (trCountArr == null) {
						trCountArr = new int[] { 0, 0, 0, 0, 0, 0 };
						trCounts.put(split[1], trCountArr);
					}
				}
				trCountArr[index] += Integer.parseInt(split[2]);
			}
			br.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void writeFeatures() {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("D:/Dennis/Uni/GoBi/A4/fdat.txt")));
			for (Entry<String, int[]> e : geneCounts.entrySet()) {
				bw.write(e.getKey() + "\n");
			}
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void writeTranscriptFeatures() {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("D:/Dennis/Uni/GoBi/A4/tr_fdat.txt")));
			for (Entry<String, TreeMap<String, int[]>> entry : transcriptCounts.entrySet()) {
				for (Entry<String, int[]> e : entry.getValue().entrySet()) {
					bw.write(e.getKey() + "\t" + entry.getKey() + "\n");
				}
			}
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void writeExpression() {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("D:/Dennis/Uni/GoBi/A4/exprs.txt")));
			for (Entry<String, int[]> e : geneCounts.entrySet()) {
				bw.write(e.getValue()[0] + "\t" + e.getValue()[1] + "\t" + e.getValue()[2] + "\t" + e.getValue()[3]
						+ "\t" + e.getValue()[4] + "\t" + e.getValue()[5] + "\n");
			}
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	public void writeTranscriptExpression() {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("D:/Dennis/Uni/GoBi/A4/tr_exprs.txt")));
			for (Entry<String, TreeMap<String, int[]>> entry : transcriptCounts.entrySet()) {
				for (Entry<String, int[]> e : entry.getValue().entrySet()) {
					bw.write(e.getValue()[0] + "\t" + e.getValue()[1] + "\t" + e.getValue()[2] + "\t" + e.getValue()[3]
							+ "\t" + e.getValue()[4] + "\t" + e.getValue()[5] + "\n");
				}
			}
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

}
