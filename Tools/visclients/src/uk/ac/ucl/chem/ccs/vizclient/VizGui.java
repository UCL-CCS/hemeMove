package uk.ac.ucl.chem.ccs.vizclient;

import info.clearthought.layout.TableLayout;

import java.awt.Container;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.MouseMotionListener;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.text.DecimalFormat;
import java.util.concurrent.ConcurrentLinkedQueue;

import javax.media.opengl.GL;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLCapabilities;
import javax.media.opengl.GLEventListener;
import javax.media.opengl.GLJPanel;
import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.BevelBorder;

import com.sun.opengl.util.Animator;

public class VizGui extends javax.swing.JPanel implements GLEventListener{
	private JTextArea notificationArea;
	static private JScrollPane jScrollPane1;
	//private GLCanvas canvas1;
	private GLJPanel canvas1;

	
	private GLCapabilities cap;

	private String hostname;
	private int port, window;
	private NetThread thread=null; 
	private TimerThread timeFramesPerSec;
	private RotateThread rotateModel;
	private Animator animator = null;
	private boolean connected = false;
	private ConcurrentLinkedQueue queue;
	private SteeringConnection nr;
	private InfoPanel ifp = null;	
	private SteeringData sd;
	private VizFrameData currentFrame = null;

	private JLabel toolTip;
	
	int pixels_x = 512;
	int pixels_y = 512;

	private long totalFramesRec = 0;
	private long totalDataRec = 0;	

	private int viewWidth;
	private int viewHeight;
	
	private double panel_width;
	private double panel_height;

	//define the views
	static public final int VIEW1 =0;
	static public final int VIEW2 =1;
	static public final int VIEW3 =2;
	static public final int VIEW4 =3;
	static public final int VIEWALL =4;
	private int view;

	private ParamWindow1 paramWindow1;
	private ParamWindow2 paramWindow2;
	private ParamWindow3 paramWindow3;
	
	Point startGlobal;
	Point endGlobal;
	
	// Default states
	boolean mousePressed = false;
	boolean mouseReleased = true;

	public VizGui(int port, String hostname, int window, Container parent) {
		super();
		queue = new ConcurrentLinkedQueue();
		this.port = port;
		this.hostname = hostname;
		this.window = window;
		initGUI();
		this.setFocusable( true );
		
		nr = new DirectBiConnection (port, hostname, window);
		sd = new SteeringData();
		view = VIEW1;
		paramWindow1 = new ParamWindow1(parent);
		paramWindow1.updateButton.addActionListener(new ActionListener () {
			public void actionPerformed(ActionEvent evt) {
				try {
					sd.setZoom_factor(Float.parseFloat(paramWindow1.zoomField.getText()));
					sd.setLatitude(Float.parseFloat(paramWindow1.latitudeField.getText()));
					sd.setLongitude(Float.parseFloat(paramWindow1.longtitudeField.getText()));
					sd.setCtr_z(Float.parseFloat(paramWindow1.vis_ctr_zField.getText()));
					sd.setCtr_y(Float.parseFloat(paramWindow1.vis_ctr_yField.getText()));
					sd.setCtr_x(Float.parseFloat(paramWindow1.vis_ctr_xField.getText()));
					send();
				} catch (Exception e) {
					System.err.println("Can parse parameters");
				}

			}
		});

		paramWindow2 = new ParamWindow2(this);
		paramWindow2.updateButton.addActionListener(new ActionListener () {
			public void actionPerformed(ActionEvent evt) {
				try {
					pixels_x =Integer.parseInt(paramWindow2.pixels_xField.getText());
					pixels_y = Integer.parseInt(paramWindow2.pixels_yField.getText());
					sd.setPixels_x(pixels_x);
					sd.setPixels_y(pixels_y);
					send();

				} catch (Exception e) {
					System.err.println("Can parse parameters");
				}
			}
		});

		paramWindow3 = new ParamWindow3(this);
		paramWindow3.updateButton.addActionListener(new ActionListener () {
			public void actionPerformed(ActionEvent evt) {
				try {
					sd.setVis_brightness(Float.parseFloat(paramWindow3.viz_brightnessField.getText()));
					sd.setVelocity_max(Float.parseFloat(paramWindow3.velocity_maxField.getText()));
					sd.setPressure_min(Float.parseFloat(paramWindow3.pressure_minField.getText()));
					sd.setPressure_max(Float.parseFloat(paramWindow3.pressure_maxField.getText()));
					sd.setStress_max(Float.parseFloat(paramWindow3.stress_maxField.getText()));
					sd.setVis_glyph_length(Float.parseFloat(paramWindow3.vis_glyph_lengthField.getText()));
					sd.setVis_streaklines_per_pulsatile_period(Float.parseFloat(paramWindow3.vis_streaklines_per_pulsatile_periodField.getText()));
					sd.setVis_streakline_length(Float.parseFloat(paramWindow3.vis_streakline_lengthField.getText()));
					send();
				} catch (Exception e) {
					System.err.println("Can parse parameters");
				}
			}
		});

	}

	public VizGui () {
		this(1, "foohost", 1024*1024, null);
	}

	public void setHostPort (int port, String hostname) {
		this.port = port;
		this.hostname = hostname;
		notificationArea.append("Changed host " + hostname + ":" + port + "\n");
		notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
	}

	public void showParamWindow1 (boolean show) {
		paramWindow1.setVisible(show);
	}

	public boolean paramWindow1Visible () {
		return paramWindow1.isVisible();
	}

	public void showParamWindow2 (boolean show) {
		paramWindow2.setVisible(show);
	}

	public boolean paramWindow2Visible () {
		return paramWindow2.isVisible();
	}

	public void showParamWindow3 (boolean show) {
		paramWindow3.setVisible(show);
	}

	public boolean paramWindow3Visible () {
		return paramWindow3.isVisible();
	}

	public void appendNotification (String message) {
		notificationArea.append(message);
		notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
	}

	public boolean startReceive() {
		if (!connected) {
			nr = new DirectBiConnection (port, hostname, window);
			connected = nr.connect();
			if (!connected) {
				notificationArea.append("Connection error: Couldn't connect to host " + hostname + ":" + port +"\n");
				notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
				return connected;
			}
			notificationArea.append("Connection started to host " + hostname + ":" + port +"\n");
			notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
			thread = new NetThread();
			thread.start();

			//start the timing thread
			timeFramesPerSec = new TimerThread();
			timeFramesPerSec.start();


			//start the animator thread
			animator = new Animator(canvas1);
			animator.setRunAsFastAsPossible(true);  //??
			animator.start();

			send();
		}
		return connected;
	}

	public boolean stopReceive() {
		if (connected){
			animator.stop();
			connected = !nr.disconnect();
			resetSteering();
			//thread.stopThread();
			timeFramesPerSec.stopThread();
			notificationArea.append("Connection terminated to host " + hostname + ":" + port +"\n");
			notificationArea.setCaretPosition(notificationArea.getDocument().getLength());
		}
		return !connected;
	}


	/** 
	 * Executed exactly once to initialize the 
	 * associated GLDrawable
	 */ 

	public void init(GLAutoDrawable drawable) {

		// print every openGL call for debugging purposes
		//drawable.setGL(new TraceGL(drawable.getGL(), System.err));
		GL gl = drawable.getGL();

		gl.glDisable (GL.GL_DEPTH_TEST);
		gl.glDisable (GL.GL_BLEND);
		gl.glShadeModel (GL.GL_FLAT);
		gl.glDisable (GL.GL_DITHER);

		gl.glClearColor( 1.0f, 1.0f, 1.0f, 1.0f ); //white
	}

	/** 
	 * Executed if the associated GLDrawable is resized
	 */ 
	public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height) {
		this.viewWidth = width;
		this.viewHeight = height;
		
		GL gl = drawable.getGL();

		gl.glViewport(0, 0, width, height);
		gl.glMatrixMode(GL.GL_PROJECTION); 
		gl.glLoadIdentity();

		gl.glOrtho(0, width, 0, height, -1, 1);
		gl.glClear(GL.GL_COLOR_BUFFER_BIT);
	}

	/** This method handles the painting of the GLDrawable */

	public void display(GLAutoDrawable drawable) {
		GL gl = drawable.getGL();
		gl.glClear(GL.GL_COLOR_BUFFER_BIT);

		VizFrameData vfd = (VizFrameData) queue.poll();
		
		if(vfd == null) {
			if(currentFrame == null) {
				// just in case - might fall through here on the first few loops.
				// won't happen if we've received at least one frame though.
				return;
			}
			// no new frame, use the current one.
			vfd = currentFrame;
		}
		else {
			// got a new frame, set it to be the current one from now on.
			currentFrame = vfd;
		}
		
		//Calendar cal = Calendar.getInstance();
		//long start_time = cal.getTimeInMillis();

		int imageWidth = vfd.getWidth();
		int imageHeight = vfd.getHeight();

		// calculate zoom required to fill panel
		float widthZoom = (float) viewWidth / (float) imageWidth;
		float heightZoom = (float) viewHeight / (float) imageHeight;
		float zoom = widthZoom < heightZoom ? widthZoom : heightZoom;

		if (view < VIEWALL) {
			gl.glPixelZoom(zoom, zoom);
			gl.glRasterPos2i((int) ((viewWidth - (imageWidth * zoom)) * 0.5f), (int) ((viewHeight - (imageHeight * zoom)) * 0.5f));
			gl.glDrawPixels(imageWidth, imageHeight, GL.GL_RGBA, GL.GL_UNSIGNED_INT_8_8_8_8, vfd.getBuffer(view));
		}
		else {
			// for viewing all frames zoom needs to be halved
			float zoomAll = zoom * 0.5f;
			gl.glPixelZoom(zoomAll, zoomAll);

			gl.glRasterPos2i((int) ((viewWidth - (imageWidth * zoom)) * 0.5f), (int) (viewHeight * 0.5));
			gl.glDrawPixels(imageWidth, imageHeight, GL.GL_RGBA, GL.GL_UNSIGNED_INT_8_8_8_8, vfd.getBuffer(0));

			gl.glRasterPos2i((int) (viewWidth * 0.5f), (int) (viewHeight * 0.5));
			gl.glDrawPixels(imageWidth, imageHeight, GL.GL_RGBA, GL.GL_UNSIGNED_INT_8_8_8_8, vfd.getBuffer(1));

			gl.glRasterPos2i((int) ((viewWidth - (imageWidth * zoom)) * 0.5f), (int) ((viewHeight - (imageHeight * zoom)) * 0.5f));
			gl.glDrawPixels(imageWidth, imageHeight, GL.GL_RGBA, GL.GL_UNSIGNED_INT_8_8_8_8, vfd.getBuffer(2));

			gl.glRasterPos2i((int) (viewWidth * 0.5f), (int) ((viewHeight - (imageHeight * zoom)) * 0.5f));
			gl.glDrawPixels(imageWidth, imageHeight, GL.GL_RGBA, GL.GL_UNSIGNED_INT_8_8_8_8, vfd.getBuffer(3));
		}

		//Calendar cal2 = Calendar.getInstance();
		//long end_time = cal2.getTimeInMillis();

		//double total_time = (end_time - start_time) / 1000.0d;
		//System.out.println("Render time: " + total_time);

		if (vfd.getVis_mouse_pressure() > -1) {
			//canvas1.setToolTipText(vfd.getVis_mouse_pressure() + " : " + vfd.getVis_stess_pressure());
			System.err.println(vfd.getVis_mouse_pressure() + " : " + vfd.getVis_stess_pressure());
			DecimalFormat df= new DecimalFormat("#####0.###"); 

			toolTip.setText(df.format(vfd.getVis_mouse_pressure()) + " : " + df.format(vfd.getVis_stess_pressure()));
		}

		if (ifp != null) {
			ifp.updatePanel (vfd.getBufferSize(), vfd.getVis_cycle(), totalFramesRec, 
					vfd.getVis_time_step(), vfd.getVis_time(), vfd.getVis_pressure_max(), vfd.getVis_pressure_min(), 
					vfd.getVis_velocity_max(), vfd.getVis_velocity_min(), vfd.getVis_stress_max(), vfd.getVis_stress_min());

			//ifp.updatePanel(vfd.getBufferSize(), 0, totalFramesRec);
			// System.err.println(vfd.getFramePerSec());
		}
	}

	/** This method handles things if display depth changes */
	public void displayChanged(GLAutoDrawable drawable,
			boolean modeChanged,
			boolean deviceChanged){
	}

	private void initGUI() {
		try {
			TableLayout thisLayout = new TableLayout(new double[][] {
					{ TableLayout.FILL },
					{ TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
						TableLayout.FILL, TableLayout.FILL,
						TableLayout.FILL, TableLayout.FILL,
						TableLayout.FILL, TableLayout.FILL,
						TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL,
						TableLayout.FILL, TableLayout.FILL,
						TableLayout.FILL, TableLayout.FILL,
						TableLayout.FILL, TableLayout.FILL,
						TableLayout.FILL } });
			this.setLayout(thisLayout);
			thisLayout.setHGap(5);
			thisLayout.setVGap(5);
			this.setMinimumSize(new java.awt.Dimension(570, 677));
			//this.setPreferredSize(new java.awt.Dimension(570, 677));
			{
				cap = new GLCapabilities(); 
				//canvas1 = new GLCanvas(cap);
				canvas1 = new GLJPanel(cap);

				this.add(canvas1, "0, 3, 0, 19");
				canvas1.addGLEventListener(this);

				toolTip = new JLabel();
				toolTip.setVisible(false);
				canvas1.add(toolTip);
				
				this.addKeyListener(new KeyListener () {
				
					public void keyPressed(KeyEvent e) { }

					public void keyReleased(KeyEvent e)  { }
				
					public void keyTyped(KeyEvent e) {
						
						System.err.println("Key pressed");

						if (e.getKeyChar() == '1') {
							view = VIEW1;

						} else if (e.getKeyChar() == '2') {
							view = VIEW2;

						} else if (e.getKeyChar() == '3') {
							view = VIEW3;

						} else if (e.getKeyChar() == '4') {
							view = VIEW4;

						} else if (e.getKeyChar() == '5') {
							view = VIEWALL;

						}
						ifp.viewChanged();
						//System.err.println("View " + view);
					}

				});

				canvas1.addMouseWheelListener(new MouseWheelListener () {

					public void mouseWheelMoved(MouseWheelEvent e) {
						float mag = sd.getZoom_factor();	
						int notches = (e.getWheelRotation());
						mag = mag - (0.1f*notches);
						sd.setZoom_factor(mag);
						send();
						// System.err.println("Mag foactor= " + sd.getZoom_factor());
					}

				});
				
				canvas1.addMouseMotionListener(new MouseMotionListener () {
					
					public void mouseDragged(MouseEvent e) {
						
						Point start = e.getPoint();
						Point end = startGlobal;
						
						if(start != null && end != null) {
						
							//System.out.println("start x pos " + start.x + " y pos " + start.y);
							//System.out.println("end   x pos " + end.x + " y pos " + end.y);
							
							if (mousePressed) {
								
								//System.out.println("start x pos " + start.x + " y pos " + start.y);

								double dx = end.getX() - start.getX();
								double dy = start.getY() - end.getY();

								panel_width = canvas1.getWidth();
								panel_height = canvas1.getHeight();
								
								dx = dx / panel_width;
								dy = dy / panel_height;
								
								dx = dx / (Math.PI / 180.0);
								dy = dy / (Math.PI / 180.0);
							
								sd.updateLongitudeLatitude(-dx, dy);
								
//								if (dx != 0) {
//									//scale 
//									dx = (90.d/panel_width)*dx;
//									sd.updateLongitude((float)dx);
//									System.err.println("dx = " + dx + " panel " + panel_width);
//									
//									
//								} if (dy != 0) {
//									dy = -1*(90.d/panel_height)*dy;
//									sd.updateLatitude((float)dy);
//									System.err.println("dy = " + dy);
//								}
								
								send();
								
							} /* else if(start != null && e.getButton() == 3) {
								Point end = e.getPoint();
								double dx = start.getX() - end.getX();
								double dy = start.getY() - end.getY();
								panel_width = canvas1.getWidth();
								panel_height = canvas1.getHeight();

								sd.setCtr_x(100*(float)(dx/panel_width));
								sd.setCtr_y(100*(float)(dy/panel_height));

								//System.err.println("Ctr x = " + sd.getCtr_x() + " crt y = " + sd.getCtr_y());
								send();
							}	 */						
	
						}
						
					}
					
					public void mouseMoved(MouseEvent e) {
						// leave blank
					}
					
				});
				
				canvas1.addMouseListener(new MouseListener () {
					Point start = null;

					public void mouseReleased(MouseEvent e) {
						
						mousePressed = false;
						mouseReleased = true;
					
						endGlobal = e.getPoint();
						
						//System.err.println("UNPRESSED");
						
						//System.err.println("Button " + e.getButton());
						
						/*if (start != null && e.getButton() == 1) {
							Point end = e.getPoint();
							double dx = start.getX() - end.getX();
							double dy = start.getY() - end.getY();

							//dddd
							panel_width = canvas1.getWidth();
							panel_height = canvas1.getHeight();

							if (dx != 0) {
								//scale 
								dx = (90.d/panel_width)*dx;
								sd.updateLongitude((float)dx);
								//System.err.println("dx = " + dx + " panel " + panel_width);
							} if (dy != 0) {
								dy = -1*(90.d/panel_height)*dy;
								sd.updateLatitude((float)dy);
								//System.err.println("dy = " + dy);
							}
							send();
						} else if(start != null && e.getButton() == 3) {
							Point end = e.getPoint();
							double dx = start.getX() - end.getX();
							double dy = start.getY() - end.getY();
							panel_width = canvas1.getWidth();
							panel_height = canvas1.getHeight();

							sd.setCtr_x(100*(float)(dx/panel_width));
							sd.setCtr_y(100*(float)(dy/panel_height));

							//System.err.println("Ctr x = " + sd.getCtr_x() + " crt y = " + sd.getCtr_y());
							send();
						} */
					}

					public void mouseClicked(MouseEvent e) {
						if (e.getButton() == 2 && view != VIEWALL) {
							Point point = e.getPoint();
							int i = (int)point.getX();	
							int j = (int)point.getY();
							
							int width = (int)canvas1.getWidth();
							int height = (int)canvas1.getHeight();

							int offset_x = (width - pixels_x)/2;
							int offset_y = (height - pixels_y)/2;
							
							j = height - j;
							
							//i = Math.round(((float)i/(float)width) * pixels_x);
							//j = Math.round(((float)j/(float)height) * pixels_y);
							
							i = i - offset_x;
							j = j - offset_y; 
						
							System.err.println("i = " + i + " j = " + j);
							
							toolTip.setLocation(point);
							toolTip.setVisible(true);					

							if (i >= 0 && i <= pixels_x && j >= 0 && j <= pixels_y) {
								sd.setVis_mouse_x(i);
								sd.setVis_mouse_y(j);
								send();
							}
							
							//System.err.println("i=" + i + " j=" + j + " x offset=" + offsetx + " y offset=" +offsety);

						}

					}

					public void mousePressed(MouseEvent e) {
						
						mousePressed = true;
						mouseReleased = false;
						
						startGlobal = e.getPoint();
						
						//System.out.println("PRESSED");
						//System.err.println("startx = " + start.getX() + " starty = " + start.getY());
					}
	
					public void mouseExited(MouseEvent e) {
					}

					public void mouseEntered(MouseEvent e) {
					}
					
				});
			}	
			{
				jScrollPane1 = new JScrollPane();
				this.add(jScrollPane1, "0, 0, 0, 2");
				jScrollPane1.setAutoscrolls(true);
				{
					notificationArea = new JTextArea();
					jScrollPane1.setViewportView(notificationArea);
					notificationArea.setText("-------Messages-------\n");
					notificationArea.setVisible(true);
					notificationArea.setBorder(BorderFactory
							.createBevelBorder(BevelBorder.LOWERED));
					notificationArea.setDoubleBuffered(true);
					notificationArea.setEditable(false);
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void viewChanged (int v) {
		view = v;
	}
	
	public int getView () {
		return view;
	}
	
	public void rotate (boolean rot) {
		if (rot) {
			//System.out.println("rotate");
			rotateModel = new RotateThread();
			rotateModel.start();
		} else {
			rotateModel.stopThread();
		}
	}

	public void kill () {
		sd.setKill(1);
		send();
	}

	public void resetSteering () {
		sd = new SteeringData();
		pixels_x = sd.getPixels_x();
		pixels_y = sd.getPixels_y();
		send();
		//	ifp.updateSteeredParams();
		view = VIEW1;
	}

	private void send() {
		if( nr == null) System.err.println("nr is NULL");
		if (nr != null && nr.isConnected()) {
			System.err.println ("Calling private void send() with (nr != null && nr.isConnected())");
			nr.send(sd);
			ifp.viewChanged();
			paramWindow1.update(sd);
			paramWindow2.update(sd);
			paramWindow3.update(sd);
			
			
			//reset some values
			if (sd.getVis_mouse_x() >= 0 || sd.getVis_mouse_y() >= 0) {
				sd.setVis_mouse_x(-1);
				sd.setVis_mouse_y(-1);
			}
			
			if (sd.getCtr_x() != 0.0f || sd.getCtr_y() != 0.0f || sd.getCtr_z() != 0.0f) {
				sd.setCtr_x(0.0f);
				sd.setCtr_y(0.0f);
				sd.setCtr_z(0.0f);

			}
			
			
		}


		//ifp.updateSteeredParams();

	}


	public void saveParameters (String filename) {
		try {
			FileOutputStream fout = new FileOutputStream(filename);
			ObjectOutputStream oos = new ObjectOutputStream(fout);
			oos.writeObject(sd);
			oos.close();
		}
		catch (Exception e) { 
			e.printStackTrace(); 
		}

	}

	public void loadParameters (String filename) {
		try {
			FileInputStream fin = new FileInputStream(filename);
			ObjectInputStream ois = new ObjectInputStream(fin);
			sd = (SteeringData) ois.readObject();
			ois.close();
			send();
		}
		catch (Exception e) { e.printStackTrace(); }

	}

	public void changeVisMode (int renView) {
		sd.setVis_mode(renView);
		send();
	}

	public int getVisMode () {
		return sd.getVis_mode();
	}

	public boolean isConnected() {
		return nr.isConnected();
	}

	public JPanel getInfoPanel () {
		ifp = new InfoPanel();
		return ifp;
	}



	//thread to rotate the model
	private class RotateThread extends Thread {
		boolean run;
		public RotateThread () {
			run=true;
		}


		public void stopThread () {
			//private boolean shouldirun = true;
			run=false;
		} 	

		public void run() {
			while(nr.isConnected() && run==true) {
				try {
					Thread.sleep(250);    			
					sd.updateLongitude(0.5f);
					send();

				} catch (Exception e) {

				}
			}


		}
	}

	//thread to time fps etc
	private class TimerThread extends Thread {
		private boolean run;
		private int toolTipDisplayTime;
		private long runTime;
		public TimerThread () {
			run=true;
			runTime = 0;
			toolTipDisplayTime = 0;
		}

		public void stopThread () {
			//private boolean shouldirun = true;
			run=false;
			runTime = 0;
			toolTipDisplayTime = 0;
		} 	

		public void run() {

			long sleepTime = 1; // seconds
			
			while(nr.isConnected() && run) {
				try {
					
					VizFrameData vfd = (VizFrameData) queue.poll();
					
					long instantDataInitial = totalDataRec;
					long instantFramesInitial = totalFramesRec;
					long instantStepInitial = vfd.getVis_time_step();
					
					long startTime = System.nanoTime();
					
					Thread.sleep(sleepTime * 1000);
					
					vfd = (VizFrameData) queue.poll();

					long instantDataFinal = totalDataRec;
					long instantFramesFinal = totalFramesRec;
					long instantStepFinal = vfd.getVis_time_step();
					
					double elapsedTime = ((System.nanoTime() - startTime) * 1.0e-9); // seconds

					double stepRate = -1.0f;
					
					if(instantStepFinal >= instantStepInitial) {
						stepRate = ((instantStepFinal-instantStepInitial) * 1.0f)/(sleepTime*1.0f);
					}
					
					runTime++;
					
					//double rate = (totalDataRec * 1.0f)/(runTime*1024.0f);
					//double frames = (totalFramesRec * 1.0f)/(runTime*1.0f);
					
					double rate = ((instantDataFinal-instantDataInitial) * 1.0f)/(elapsedTime*1024.0f);
					double frames = ((instantFramesFinal-instantFramesInitial) * 1.0f)/(elapsedTime*1.0f);
					
					System.out.println("Frame rate " + frames + " fps, Data rate " + rate + " KB/s " + " Step rate " + stepRate + " steps/s");					
					
					ifp.updateFPS(frames, rate);

					//System.err.println("Total data " + totalDataRec + " Total frames " + totalFramesRec + " Time " + runTime);
				/*	if (toolTip.isVisible() && toolTipDisplayTime == 0) {
						toolTipDisplayTime = 5;
					} else if (toolTipDisplayTime > 0) {
						toolTipDisplayTime--;
						if (toolTipDisplayTime == 0) {
							toolTip.setText("");
							toolTip.setVisible(false);
						}
					}*/
					
				} catch (Exception e) {

				}
			}

		}

	}

	//thread to do the network receive
	private class NetThread extends Thread {

		//TimerThread timeFramesPerSec;

		public NetThread () {

		}

		public void stopThread () {
			nr.disconnect();
			//timeFramesPerSec.stop();
		}

		public void run() {
			//timeFramesPerSec.start();

			while(nr.isConnected()) {
				VizFrameData vfd = nr.getFrame();
				if (vfd != null) {
					//System.out.println("Frames stacked " + queue.size());
					if( queue.size() > 30 ) {
						// Removed the oldest frame from the queue
						queue.poll();
						//System.out.println("Removing frame..." + queue.size());						
					}
					queue.offer(vfd);
					totalFramesRec++;
					totalDataRec = totalDataRec + vfd.getBufferSize();
				}
			}
			stopReceive();

		}
	}

	//side info panel 
	private class InfoPanel extends JPanel {

		private JLabel jLabel4;
		private JTextField dataRate;
		private JLabel jLabel5;
		private JLabel jLabel6;
		private JTextField timeStep;
		private JTextField viewing;
		private JTextField frameSize;
		private JTextField frameRecNo;
		private JTextField cycle;
		private JLabel jLabel3;
		private JLabel jLabel2;
		private JLabel jLabel1;
		private JTextField framePerSec;
		private JLabel jLabel10;
		private JLabel jLabel7;
		private JTextField visTime;
		
		private JLabel jLabel8;
		private JTextField pressure;
		private JLabel jLabel9;
		private JTextField vel;
		private JLabel jLabel11;
		private JTextField stress;
		
		public InfoPanel () {
			super();

			TableLayout infoPanelLayout = new TableLayout(new double[][] {{TableLayout.FILL}, {TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL, TableLayout.FILL}});
			infoPanelLayout.setHGap(5);
			infoPanelLayout.setVGap(5);
			this.setLayout(infoPanelLayout);
			this.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
			{
				jLabel2 = new JLabel();
				this.add(jLabel2, "0, 0");
				jLabel2.setText("pulsatile cycle id");
				jLabel2.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				cycle = new JTextField();
				this.add(cycle, "0, 1");
			}
			{
				jLabel1 = new JLabel();
				this.add(jLabel1, "0, 2");
				jLabel1.setText("# frames received");
				jLabel1.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
				jLabel1.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				frameRecNo = new JTextField();
				this.add(frameRecNo, "0, 3");
				frameRecNo.setText("");
			}
			{
				jLabel3 = new JLabel();
				this.add(jLabel3, "0, 4");
				jLabel3.setText("Frame size (bytes)");
				jLabel3.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				frameSize = new JTextField();
				this.add(frameSize, "0, 5");
				frameSize.setText("");
			}
			{
				jLabel5 = new JLabel();
				this.add(jLabel5, "0, 6");
				jLabel5.setText("HemeLB time step");
				jLabel5.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				timeStep = new JTextField();
				this.add(timeStep, "0, 7");
				timeStep.setText("");
			}
			{
				jLabel4 = new JLabel();
				this.add(jLabel4, "0, 8");
				jLabel4.setText("Data rate (KB/s)");
				jLabel4.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				dataRate = new JTextField();
				this.add(dataRate, "0, 9");
			}
			{
				jLabel6 = new JLabel();
				this.add(jLabel6, "0, 10");
				jLabel6.setText("View");
				jLabel6.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				viewing = new JTextField();
				this.add(viewing, "0, 11");			
				viewing.setText("");
			}
			{
				jLabel10 = new JLabel();
				this.add(jLabel10, "0, 12");
				jLabel10.setText("fps");
				jLabel10.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				framePerSec = new JTextField();
				this.add(framePerSec, "0, 13");
			}
			{
				jLabel7 = new JLabel();
				this.add(jLabel7, "0, 14");
				jLabel7.setText("Intra cycle time (s)");
				jLabel7.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				visTime = new JTextField();
				this.add(visTime, "0, 15");
			}

			{
				jLabel8 = new JLabel();
				this.add(jLabel8, "0, 16");
				jLabel8.setText("Pressure min/max (mm.Hg)");
				jLabel8.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				pressure = new JTextField();
				this.add(pressure, "0, 17");
			}
			{
				jLabel9 = new JLabel();
				this.add(jLabel9, "0, 18");
				jLabel9.setText("Velocity min/max (m/s)");
				jLabel9.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				vel = new JTextField();
				this.add(vel, "0, 19");
			}
			{
				jLabel11 = new JLabel();
				this.add(jLabel11, "0, 20");
				jLabel11.setText("Stress min/max (Pa)");
				jLabel11.setFont(new java.awt.Font("Dialog",1,12));
			}
			{
				stress = new JTextField();
				this.add(stress, "0, 21");
			}
			
		}


		public void updateFPS (double fps, double bps) {
			DecimalFormat df= new DecimalFormat("#####0.###"); 
			framePerSec.setText(df.format(fps));
			dataRate.setText(df.format(bps));
		}

		public void viewChanged () {
			switch (view) {
			case 0:
				viewing.setText("Velocity Volume Rendering");
				break;	
			case 1:
				viewing.setText("Stress Volume Rendering");
				break;
			case 2:
				viewing.setText("Wall Pressure");
				break;
			case 3:
				viewing.setText("Wall Stress");
				break;
			case 4:
				viewing.setText("All Views");
				break;
			}

		}

		public void updatePanel (int fSize, int cyc, long fRecNo, int timeS, double visT, double pMax, double pMin, double vMax, double vMin, double sMax, double sMin) {
			
			cycle.setText(Integer.toString(cyc));
			frameRecNo.setText(Long.toString(fRecNo));
			frameSize.setText(Integer.toString(fSize));
			timeStep.setText(Integer.toString(timeS));
			
			DecimalFormat df2 = new DecimalFormat("#####0.####");
			visTime.setText(df2.format(visT));
			
			DecimalFormat df = new DecimalFormat("#####0.###");
			pressure.setText(df.format(pMin) + "/" + df.format(pMax));
			vel.setText(df.format(vMin) + "/" + df.format(vMax));
			stress.setText(df.format(sMin) + "/" + df.format(sMax));

		}

	}

}
