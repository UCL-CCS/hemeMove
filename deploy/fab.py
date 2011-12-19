"""
Fabric definitions for HemeLB
Usage:
	With current working directory anywhere inside the HemeLB mercurial checkout, execute
		"fab <machinename> <task>"
	for example:
		"fab hector deploy_cold"
	
	Before use, you MUST copy deploy/machines_user_example.json as deploy/machines_user.json and fill in your personal details.
	For smoothest usage, you should also ensure ssh keys are in place so the target machine can access the mercurial repository
	of hemelb.
"""

from templates import *
from machines import *
from fabric.contrib.project import *
import time

@task
def clone():
	"""Checkout the repository afresh on the remote machine, wiping out any existing checkout."""
	run(template("mkdir -p $remote_path"))
	with cd(env.remote_path):
		run(template("rm -rf $repository"))
		run(template("hg clone $hg/$repository"))

@task(alias='cold')
def deploy_cold():
	"""Completely checkout, build, and install hemelb, from new, on the remote machine."""
	execute(clone)
	execute(clear_build)
	execute(prepare_paths)
	execute(configure)
	execute(build)
	execute(install)
	execute(build_python_tools)

@task
def update_build():
	"""Update the remote repository and perform an incremnental build and install"""
	execute(update)
	execute(require_recopy)
	execute(configure)
	execute(build)
	execute(install)
	execute(build_python_tools)

@task 
def require_recopy():
	"""Notify the build system that the HemeLB C++ code has changed."""
	run(template("touch $build_path/hemelb-prefix/src/hemelb-stamp/hemelb-mkdir"))

@task
def update():
	"""Update the remote mercurial repository"""
	with cd(env.repository_path):
		run("hg pull")
		run("hg update")

@task 
def prepare_paths():
	"""Prepare remote locations to store scripts, results, and config files"""
	run(template("mkdir -p $scripts_path"))
	run(template("mkdir -p $results_path"))
	run(template("mkdir -p $config_path"))
	
@task
def clear_build():
	"""Wipe out existing built and installed HemeLB."""
	run(template("rm -rf $build_path"))
	run(template("rm -rf $install_path"))
	run(template("mkdir -p $build_path"))
	run(template("mkdir -p $install_path"))
	

@task
def clean():
	"""Clean remote build area."""
	with cd(env.build_path):
		run("make clean")

@task(alias='tools')
def build_python_tools():
	"""Build and install python scripts used for working with HemeLB"""
	with cd(env.tools_path):
		run("python setup.py build")

@task
def stat():
	"""Check the remote message queue status"""
	#TODO: Respect varying remote machine queue systems.
	run("qstat")
	
@task
def monitor():
	"""Continuously check and report on the remote queue status, ctrl-C to interrupt"""
	while True:
		time.sleep(30)
		execute(stat)

@task
def configure():
	"""Execute the CMake configure step for HemeLB and dependencies."""
	with cd(env.build_path):
		with prefix(env.build_prefix):
			run(template(
			"cmake $repository_path -DCMAKE_INSTALL_PREFIX=$install_path "+
			"-DDEPENDENCIES_INSTALL_PATH=$install_path $cmake_flags"
			))

@task
def build(verbose=False):
	"""Execute the CMake build step for HemeLB and dependencies."""
	with cd(env.build_path):
		run(template("rm -rf hemelb_prefix/build"))
		with prefix(env.build_prefix):
			if verbose:
				run("make VERBOSE=1")
			else:
				run("make")

@task
def install():
	"""Execute the CMake install step for HemeLB and dependencies."""
	with cd(env.build_path):
		with prefix(env.build_prefix):
			run("make install")
			run(template("chmod u+x $install_path/bin/unittests_hemelb $install_path/bin/hemelb"))
		
@task
def revert(args="--all"):
	"""Revert local changes in the remote repository.
	 Including those made through 'fab ... patch' and 'fab ... sync'.
	Specify a path relative to the repository root to revert only some files or directories
	"""
	with cd(env.repository_path):
		run("hg revert %s"%args)
	
@task	
def send_distributions():
	"""Transmit dependency tarballs from local dependencies/distributions to remote.
	Useful to prepare a build on target machines with CMake before 2.8.4, where
	HTTP redirects are not followed.
	"""
	put(os.path.join(env.localroot,'dependencies','distributions','*.tar.gz'),
		env.pather.join(env.repository_path,'dependencies','distributions'))

@task
def sync():
	"""Update the remote repository with all changes from the local repository, via rsync.
	Respects the local .hgignore files to avoid sending unnecessary information.
	"""
	rsync_project(
		remote_dir=env.repository_path,
		local_dir=env.localroot+'/',
		exclude=map(lambda x: x.replace('\n',''),
		list(open(os.path.join(env.localroot,'.hgignore')))+
		['.hg']+
		list(open(os.path.join(env.localroot,'RegressionTests','.hgignore')))
		)
	)

@task
def patch(args=""):
	"""Update the remote repository with all changes from the local repository, via hg diff.
	Specify a path relative to the repository root to patch only some files or directories
	e.g 'fab legion patch Code/main.cc'
	"""
	local("hg diff %s> fabric.diff"%args)
	put("fabric.diff",env.pather.join(env.remote_path,env.repository))
	with cd(env.repository_path):
		run("patch -p1 < fabric.diff")

def with_job(name):
	"""Augment the fabric environment with information regarding a particular job name.
	Definitions created:
	job_results: the remote location where job results should be stored
	job_results_local: the local location where job results should be stored
	"""
	env.job_results=env.pather.join(env.results_path,name)
	env.job_results_local=os.path.join(env.local_results,name)
	env.job_results_contents=env.pather.join(env.job_results,'*')
	env.job_results_contents_local=os.path.join(env.job_results_local,'*')


def with_config(name):
	"""Internal: augment the fabric environment with information regarding a particular configuration name.
	Definitions created:
	job_config_path: the remote location where the config files for the job should be stored
	job_config_path_local: the local location where the config files for the job may be found
	"""
	env.job_config_path=env.pather.join(env.config_path,name)
	env.job_config_path_local=os.path.join(env.local_configs,name)
	env.job_config_contents=env.pather.join(env.job_config_path,'*')
	env.job_config_contents_local=os.path.join(env.job_config_path_local,'*')
	

@task
def fetch_configs(config=''):
	"""
	Fetch config files from the remote, via rsync.
	Specify a config directory, such as 'cylinder' to copy just one config.
	Config files are stored as, e.g. cylinder/config.dat and cylinder/config.xml
	Local path to use is specified in machines_user.json, and should normally point to a mount on entropy,
	i.e. /store4/blood/username/config_files
	This method is not intended for normal use, but is useful when the local machine cannot have an entropy mount,
	so that files can be copied to a local machine from entropy, and then transferred to the compute machine,
	via 'fab entropy fetch_configs; fab legion put_configs'
	"""
	with_config(config)
	local(template("rsync -pthrvz $username@$remote:$job_config_path/ $job_config_path_local"))

@task
def put_configs(config=''):
	"""
	Transfer config files to the remote, for use in launching jobs, via rsync.
	Specify a config directory, such as 'cylinder' to copy just one configuration.
	Config files are stored as, e.g. cylinder/config.dat and cylinder/config.xml
	Local path to find config directories is specified in machines_user.json, and should normally point to a mount on entropy,
	i.e. /store4/blood/username/config_files
	If you can't mount entropy, 'fetch_configs' can be useful, via 'fab entropy fetch_configs; fab legion put_configs'
	"""
	with_config(config)
	run(template("mkdir -p $job_config_path"))
	rsync_project(local_dir=env.job_config_path_local+'/',remote_dir=env.job_config_path)

@task
def put_results(name=''):
	"""
	Transfer result files to a remote.
	Local path to find result directories is specified in machines_user.json.
	This method is not intended for normal use, but is useful when the local machine cannot have an entropy mount,
	so that results from a local machine can be sent to entropy, via 'fab legion fetch_results; fab entropy put_results'
	"""
	with_job(name)
	run(template("mkdir -p $job_results"))
	rsync_project(local_dir=env.job_results_local+'/',remote_dir=env.job_results)
	
@task
def fetch_results(name=''):
	"""
	Fetch results of remote jobs to local result store.
	Specify a job name to transfer just one job.
	Local path to store results is specified in machines_user.json, and should normally point to a mount on entropy,
	i.e. /store4/blood/username/results.
	If you can't mount entropy, 'put results' can be useful,  via 'fab legion fetch_results; fab entropy put_results'
	"""
	with_job(name)
	local(template("rsync -pthrvz $username@$remote:$job_results/ $job_results_local"))

@task
def clear_results(name=''):
	"""Completely clear and wipe all result files from the remote."""
	with_job(name)
	run(template('rm -rf $job_results_contents'))		

@task
def test():
	"""Submit a unit-testing job to the remote queue."""
	execute(job,script='unittests',name='unittests-$build_number-$machine_name',nodes=1)
		
@task
def hemelb(**args):
	"""Submit a HemeLB job to the remote queue.
	The job results will be stored with a name pattern of $config-$build_number-$machine_name-$nodes,
	e.g. cylinder-abcd1234-legion-256
	Keyword arguments:
		config : config directory to use to define geometry, e.g. config=cylinder
		nodes : number of compute nodes to request
		images : number of images to take
		snapshots : number of snapshots to take
		steering : steering session i.d.
		wall_time : wall-time job limit
		memory : memory per node
	"""
	options=dict(script='hemelb',
		name='$config-$build_number-$machine_name-$nodes',
		nodes=4,images=10, snapshots=10, steering=1111, wall_time='0:15:0',memory='2G')
	options.update(args)
	execute(put_configs,args['config'])
	execute(job,**options)

@task(alias='regress')
def regression_test():
	"""Submit a regression-testing job to the remote queue."""
	execute(job,script='regression',name='regression-$build_number-$machine_name')

@task
def job(**args):
	"""Internal: Execute a generic job on the remote machine. Use hemelb, regress, or test instead."""
	env.update(name=args['script'],wall_time='0:1:0',nodes=4,memory='1G') #defaults
	env.update(**args)
	env.update(name=template(env.name))
	with_job(env.name)
	
	script_name=template("$machine_name-$script")
	env.job_script=script_template(script_name)
	env.dest_name=env.pather.join(env.scripts_path,env.pather.basename(env.job_script))
	put(env.job_script,env.dest_name)
	
	run(template("mkdir -p $job_results"))
	run(template("cp $dest_name $job_results"))
	run(template("cp $build_cache $job_results"))
	run(template("chmod u+x $dest_name"))
	run(template("$job_dispatch $dest_name"))
	