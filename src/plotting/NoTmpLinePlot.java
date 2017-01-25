package plotting;

import java.io.File;
import java.util.Vector;

import debugStuff.DebugMessageFactory;
import io.AllroundFileWriter;
import io.ConfigReader;
import javafx.util.Pair;

public class NoTmpLinePlot extends Plot {

	Vector<Vector<Object>> x;
	Vector<Vector<Object>> y;

	Vector<Object> legendLabels;

	double maxX = 0, maxY = 0;
	double minX = 0, minY = 0;

	public boolean showLegend = false;

	public String filename;
	private String outputPath, type;

	public NoTmpLinePlot(Pair<Vector<Vector<Object>>, Vector<Vector<Object>>> pair, String title, String xLab,
			String yLab, double maxX, double maxY, boolean logScaleX, boolean logScaleY, String outputPath,
			String type) {

		this.outputPath = outputPath;
		this.type = type;

		this.maxX = maxX;
		this.maxY = maxY;

		if (logScaleX) {
			this.x = logScaleX(pair.getKey());

			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}

		} else {
			this.x = pair.getKey();
		}

		if (logScaleY) {
			this.y = logScaleY(pair.getValue());

			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}

		} else {
			this.y = pair.getValue();
		}

		setTitle(title);
		filename = title;
		setXLab(logScaleX ? xLab + " (log10)" : xLab);
		setYLab(logScaleY ? yLab + " (log10)" : yLab);

		legendLabels = new Vector<>();
		if (showLegend) {
			for (int i = 0; i < x.size(); i++) {
				legendLabels.add(String.valueOf(i));
			}
		}
	}

	public Vector<Vector<Object>> logScaleX(Vector<Vector<Object>> in) {

		Vector<Vector<Object>> tmp = new Vector<>();
		double max = 0;

		for (Vector<Object> v : in) {
			Vector<Object> tmp2 = new Vector<>();
			for (Object v2 : v) {
				double x = Math.log10((double) ((int) v2 + 1.0));
				tmp2.add(x);
				max = Math.max(max, x);
			}
			tmp.add(tmp2);
		}

		this.maxX = (int) (max + 1);
		return tmp;
	}

	public Vector<Vector<Object>> logScaleY(Vector<Vector<Object>> in) {

		Vector<Vector<Object>> tmp = new Vector<>();
		double max = 0;

		for (Vector<Object> v : in) {
			Vector<Object> tmp2 = new Vector<>();
			for (Object v2 : v) {
				double x = Math.log10((double) ((int) v2 + 1.0));
				tmp2.add(x);
				max = Math.max(max, x);
			}
			tmp.add(tmp2);
		}

		this.maxY = (int) (max + 1);
		return tmp;
	}

	public void plot() {
		RExecutor r = new RExecutor(generateCommand(
				ConfigReader.readConfig().get("output_directory") + "/" + this.outputPath + "_" + type + ".png"));

		Thread t = new Thread(r);
		t.start();

		try {
			DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE, "Wait for R to plot..");
			t.join();
			DebugMessageFactory.printNormalDebugMessage(ConfigReader.DEBUG_MODE,
					"(" + this.filename + ".png) R thread terminated.");

		} catch (InterruptedException e) {
			throw new RuntimeException("R did not exit properly!");
		}

	}

	public void addLegendVector(Vector<Object> vector) {
		showLegend = true;
		this.legendLabels = vector;
	}

	@Override
	public String generateCommand(String filename) {

		File tmp = new File(
				ConfigReader.readConfig().get("output_directory") + "/" + this.outputPath + "_" + type + ".txt");

		AllroundFileWriter.writeVector(x.get(0), tmp);
		AllroundFileWriter.writeVector(y.get(0), tmp, true);

		for (int i = 1; i < x.size(); i++) {
			AllroundFileWriter.writeVector(x.get(i), tmp, true);
			AllroundFileWriter.writeVector(y.get(i), tmp, true);
		}
		if (showLegend) {
			AllroundFileWriter.writeVector(this.legendLabels, tmp, true);
		}

		String command = "";
		command += String.format("png(\'%s\',width=3.25,height=3.25,units=\'in\',res=400,pointsize=4);", filename);
		command += String.format("x<-scan(\'%s\',nlines=1,skip=0);", tmp.getAbsolutePath().replace("\\", "/"));
		command += String.format("y<-scan(\'%s\',nlines=1,skip=1);", tmp.getAbsolutePath().replace("\\", "/"));

		if (showLegend) {
			command += String.format("par(mar=c(%s, %s, %s, %s), xpd = TRUE);", "5.1", "4.1", "4.1", "9.1");
		} else {
			command += String.format("par(mar=c(%s, %s, %s, %s), xpd = TRUE);", "5.1", "4.1", "4.1", "4.1");
		}

		command += String.format("plot(x,y,ann=F,type=\'l\',xlim=range(" + minX + ":" + maxX + "),ylim=range(" + minY
				+ ":" + maxY + "),col=1);");

		int counter = 2;

		for (int i = 1; i < x.size(); i++) {
			command += String.format("x" + i + "<-scan(\'%s\',nlines=1,skip=" + (counter) + ");",
					tmp.getAbsolutePath().replace("\\", "/"));
			counter += 1;
			command += String.format("y" + i + "<-scan(\'%s\',nlines=1,skip=" + (counter) + ");",
					tmp.getAbsolutePath().replace("\\", "/"));
			counter += 1;
			command += String.format("lines(x" + i + ",y" + i + ",type=\'l\',col=" + (i + 1) + ");",
					tmp.getAbsolutePath().replace("\\", "/"));

		}

		if (showLegend) {
			command += String.format("ll<-scan(\'%s\',nlines=1,skip=" + (counter) + ",what=character());",
					tmp.getAbsolutePath().replace("\\", "/"));
			command += String
					.format("legend(\'topright\', inset=c(-0.25, 0), legend=ll, col=1:" + x.size() + ", lty=c(1:1));");
		}
		command += String.format("title(main=\'%s\', xlab=\'%s\', ylab=\'%s\');", super.title, super.xLab, super.yLab);
		command += "dev.off();";

		return command;

	}
}
