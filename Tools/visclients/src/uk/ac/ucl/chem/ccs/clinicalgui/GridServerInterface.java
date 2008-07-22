package uk.ac.ucl.chem.ccs.clinicalgui;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Vector;

import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import org.globus.ftp.FileInfo;

/**
 * @author Konstantin Voevodski
 *
 * implements methods for putting and getting data sets from the grid ftp server;
 * grid ftp server stores the model data generated by the segmentation tool:
 * each "data set" contains a config and pars file (used as input by hemelb), 
 * as well as a note file; data sets are organized by (patient,study,series,date/time)
 */
public class GridServerInterface {
	
	private static String HOST = "bunsen.chem.ucl.ac.uk";
	private static String ROOT = "/home/konstantin/models";
	public static String topString = "Available hemelb input";
	
	public static void uploadFiles(String path1, String path2, String path3, String patientId, String studyId, String seriesId, String dateTime) throws Exception{
	//uploads the three files associated with a particular model (two output files of seg tool plus one note file) to the grid ftp server
	//models are stored on the server organized by (patient,study,series,date)
		addDirectory(ROOT, patientId);
		addDirectory(ROOT + "/" + patientId, studyId);
		addDirectory(ROOT + "/" + patientId + "/" + studyId, seriesId);
		addDirectory(ROOT + "/" + patientId + "/" + studyId + "/" + seriesId, dateTime);
		uploadToDirectory(ROOT + "/" + patientId + "/" + studyId + "/" + seriesId + "/" + dateTime, path1, "config.dat");
		uploadToDirectory(ROOT + "/" + patientId + "/" + studyId + "/" + seriesId + "/" + dateTime, path2, "pars.asc");
		uploadToDirectory(ROOT + "/" + patientId + "/" + studyId + "/" + seriesId + "/" + dateTime, path3, "notes.txt");	
	}

	
	private static void addDirectory(String workingDir, String dir) throws Exception {
	//creates dir in workingDir on grid ftp server, if it is not already there
		MyGridFtp ftp = new MyGridFtp(HOST, 2811);
		ftp.changeDir(workingDir);
		Vector<FileInfo> v = ftp.list();
        for(FileInfo f : v){
        	if(f.getName().equals(dir))
        		return;
        }
        ftp.makeDir(dir);
        ftp.close();
	}
	
	private static void uploadToDirectory(String workingDir, String localFile, String remoteFile) throws Exception{
	//uploads localFile to remoteFile in workingDir on grid ftp server
		MyGridFtp ftp = new MyGridFtp(HOST, 2811);
		ftp.changeDir(workingDir);
		ftp.upload(localFile, remoteFile);
		ftp.close();
	}
	
	private static Vector<String> listDirectory(String workingDir) throws Exception{
	//lists contents of workingDir: subroutine of populateServerData()
		Vector<String> o = new Vector<String>();
		MyGridFtp ftp = new MyGridFtp(HOST, 2811);
		ftp.changeDir(workingDir);
		Vector<FileInfo> v = ftp.list();
        for(FileInfo f : v){
		    String name = f.getName();
		    if(name.equals(".") || name.equals(".."))
		    	continue;
        	o.add(name);
        	//System.out.println(name);
        }
        ftp.close();
        return o;
	}
	
	/*
	public static DefaultMutableTreeNode populateServerData() throws Exception {
		DefaultMutableTreeNode top = new DefaultMutableTreeNode(topString);
		Vector<String> patients = listDirectory(ROOT);
		for(String patient : patients){
			DefaultMutableTreeNode patientNode = new DefaultMutableTreeNode(patient);
			top.add(patientNode); //add patient node
			Vector<String> studies = listDirectory(ROOT + "/" + patient);
			for(String study : studies){
				DefaultMutableTreeNode studyNode = new DefaultMutableTreeNode(study);
				patientNode.add(studyNode); //add study node
				Vector<String> series = listDirectory(ROOT + "/" + patient + "/" + study);
				for(String ser : series){
					DefaultMutableTreeNode seriesNode = new DefaultMutableTreeNode(ser);
					studyNode.add(seriesNode); //add series node
					Vector<String> dates = listDirectory(ROOT + "/" + patient + "/" + study + "/" + ser);
					for(String date : dates){
						DefaultMutableTreeNode dateNode = new DefaultMutableTreeNode(date);
						seriesNode.add(dateNode); //add date node
						Vector<String> files = listDirectory(ROOT + "/" + patient + "/" + study + "/" + ser + "/" + date);
						for(String file : files){
							DefaultMutableTreeNode fileNode = new DefaultMutableTreeNode(file);
							dateNode.add(fileNode); //add file node
						}
					}
				}
			}
			
		}
		return top;
	}
	*/
	
	
	public static JTree populateServerData(String[] highlightedPath) throws Exception {
	//returns a tree representing the file structure of the grid ftp server
		DefaultMutableTreeNode top = new DefaultMutableTreeNode(topString);
		DefaultMutableTreeNode[] hp = new DefaultMutableTreeNode[5];
		hp[0] = top;
		Vector<String> patients = listDirectory(ROOT);
		for(String patient : patients){
			DefaultMutableTreeNode patientNode = new DefaultMutableTreeNode(patient);
			if(highlightedPath !=  null && patient.equals(highlightedPath[0]))
				hp[1] = patientNode;
			top.add(patientNode); //add patient node
			Vector<String> studies = listDirectory(ROOT + "/" + patient);
			for(String study : studies){
				DefaultMutableTreeNode studyNode = new DefaultMutableTreeNode(study);
				if(highlightedPath != null && study.equals(highlightedPath[1]))
					hp[2] = studyNode;
				patientNode.add(studyNode); //add study node
				Vector<String> series = listDirectory(ROOT + "/" + patient + "/" + study);
				for(String ser : series){
					DefaultMutableTreeNode seriesNode = new DefaultMutableTreeNode(ser);
					if(highlightedPath != null && ser.equals(highlightedPath[2]))
						hp[3] = seriesNode;
					studyNode.add(seriesNode); //add series node
					Vector<String> dates = listDirectory(ROOT + "/" + patient + "/" + study + "/" + ser);
					for(String date : dates){
						DefaultMutableTreeNode dateNode = new DefaultMutableTreeNode(date);
						if(highlightedPath != null && date.equals(highlightedPath[3]))
							hp[4] = dateNode;
						seriesNode.add(dateNode); //add date node
						Vector<String> files = listDirectory(ROOT + "/" + patient + "/" + study + "/" + ser + "/" + date);
						for(String file : files){
							DefaultMutableTreeNode fileNode = new DefaultMutableTreeNode(file);
							dateNode.add(fileNode); //add file node
						}
					}
				}
			}
			
		}
		JTree tree = new JTree(top);
		if(highlightedPath != null)
			tree.setSelectionPath(new TreePath(hp));
		return tree;
	}
	
	public static String getModelNote(String patientId, String studyId, String seriesId, String dateTime) throws Exception{
	//downloads the note file corresponding to a particular dataset on the grid ftp server, and returns its contents as a string
		MyGridFtp ftp = new MyGridFtp(HOST, 2811);
		ftp.changeDir(ROOT + "/" + patientId + "/" + studyId + "/" + seriesId + "/" + dateTime);
		ftp.download("/tmp/notes.txt", "notes.txt");
		ftp.close();
		BufferedReader reader = new BufferedReader(new FileReader("/tmp/notes.txt"));
		return reader.readLine();
	}
}