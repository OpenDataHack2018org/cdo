require 'json'
require 'net/ssh'
require 'open3'
require 'tempfile'
require 'pp'
require 'colorize'
require 'parallel'
require 'logger'
# configuration {{{ ============================================================
# load user setting if available
RC          = "#{ENV['HOME']}/.rake.json"
@userConfig = ( File.exist?(RC) ) ? JSON.load(File.open(RC)) : {}
if @userConfig.empty? then
  warn "No host information!!"
  exit(1)
end
# get setup from the environment
@debug      = Rake.application.options.silent ? false : true
@user       = ENV['USER']
# internal variables
@_help      = {}
# }}}

# helper methods {{{ ===========================================================
# general debugging output
def dbg(msg); pp msg if @debug; end

# return name of current branch
def getBranchName; `git branch`.split("\n").grep(/^\*/)[0].split[-1].tr(')',''); end

# wrapper for calling commands locally
#
# stdout is shown in debug mode only
# stderr is always shown
def executeLocal(cmd)
  system(cmd)
  return
  Open3.popen3(cmd) {|stdin, stdout, stderr, external|
    # read from stdout and stderr in parallel
    { :out => stdout, :err => stderr }.each {|key, stream|
      Thread.new do
        until (line = stream.gets).nil? do
          puts line.chomp                       if :out == key
          puts line.chomp.colorize(color: :red) if :err == key
        end
      end
    }
    # Don't exit until the external process is done
    external.join
  }
end

# execute remote command and collect
#
# stdout is shown in debug mode only
# stderr is always shown
def executeRemote(command, builder)
  Net::SSH.start(builder.hostname,builder.username,
                 :config => true, :compression => true) do |ssh|

    stdout_data = ""
    stderr_data = ""
    exit_code   = nil
    exit_signal = nil

    ssh.open_channel do |channel|
      channel.exec(command) do |ch, success|
        unless success
          raise "FAILED: couldn't execute command #{command}"
        end

        channel.on_data do |ch, data|
          stdout_data += data
          $stdout.write(data) if @debug
        end

        channel.on_extended_data do |ch, type, data|
          stderr_data += data
          $stderr.write(data)# if @debug
        end

        channel.on_request("exit-status") do |ch, data|
          exit_code = data.read_long
        end

        channel.on_request("exit-signal") do |ch, data|
          exit_signal = data.read_long
        end
      end
    end
    ssh.loop
  end
end
#
# execution wrapper
def execute(command, builder)
  # 1)work in the target directory, ONLY
  # 2)load the user given config files for environment setup
  commands = (builder.envConfigFiles.map {|rcfile| "test -f #{rcfile} && source #{rcfile}"} +
            ["test ! -d #{builder.targetDir} && mkdir -p #{builder.targetDir}",
             "cd #{builder.targetDir}",
             command]).join(';')

  dbg(commands)

  if builder.isLocal? then
    executeLocal(commands)
  else
    executeRemote(commands,builder)
  end
end
#
# synchronization for a given host
def doSync(builder)
  # make sure, that the remote target directory is present
  execute("mkdir -p #{builder.targetDir}", builder)

  # basic rsync options
  # * remove everything on the remote site first
  rsyncOpts = "--delete-excluded --delete"
  # * keep old stuff on the remote site
  rsyncOpts = "-L"

  # collect the source files
  file = Tempfile.new("rsyncCdoTempfiles2Transfer")
  begin
    file.write(`git ls-files`)
    file.write(`git submodule foreach 'git ls-files | sed "s|^|$path/|"'`)

    # call rsync for a given host
    if builder.isLocal?
      syncCmd = "rsync #{rsyncOpts} -avz --files-from=#{file.path} . #{builder.targetDir}"
    else
      syncCmd = "rsync #{rsyncOpts} -avz --files-from=#{file.path}  -e ssh . #{builder.username}@#{builder.hostname}:#{builder.targetDir}"
    end
    dbg(syncCmd)
    executeLocal(syncCmd)
  ensure
    file.close
    file.unlink
  end
end
#
# construct task from builder object
def builder2task(builder,useHostAsName=false,syncSource=true)
  baseTaskName = useHostAsName ? builder.host : "#{builder.host}#{builder.compiler.upcase}"
  toDo         = lambda {|what| "#{baseTaskName}_#{what}".to_sym}

  if syncSource then
    @_help[:sync] = "sync source files " unless @_help.has_key?(:sync)
    task toDo[:sync] do |t|
      dbg("sync source  code for branch:" + getBranchName)
      doSync(builder)
    end
  end

  @_help[:conf]= \
    "run configure on host: ./config/default with user settings activated" unless @_help.has_key?(:conf)
  task toDo[:conf] do |t|
    dbg("call #{builder.configureCall}")
    execute("#{builder.configureCall}",builder)
  end

  @_help[:make] = \
    "run 'make'" unless @_help.has_key?(:make)
  task toDo[:make].to_sym do |t|
    execute("make -j4",builder)
  end

  @_help[:check] = \
    "run 'make check'" unless @_help.has_key?(:check)
  task toDo[:check] do |t|
    execute("make check",builder)
  end

  @_help[:checkSerial] = \
    "check with serialized IO (-L option)" unless @_help.has_key?(:checkSerial)
  task toDo[:checkSerial] do |t|
    execute("make check CDO='#{builder.targetDir}/src/cdo -L'",builder)
  end


  @_help[:clean] = \
    "run 'make clean'" unless @_help.has_key?(:clean)
  task toDo[:clean] do |t|
    execute("make clean",builder)
  end
  @_help[:cleanSync] = \
    "rm target source dir and perform a fresh sync" unless @_help.has_key?(:cleanSync)
  task toDo[:cleanSync] do |t|
    execute("cd ; rm -rf #{builder.targetDir}",builder)
    doSync(builder) if syncSource
  end

  @_help[:checkV] = \
    "run './src/cdo -V' " unless @_help.has_key?(:checkV)
  task toDo[:checkV] do |t|
    execute("./src/cdo -V",builder)
  end

  @_help[:showLog] = "show config.log file" unless @_help.has_key?(:make)
  task toDo[:showLog] do |t|
    execute("cat config.log",builder)
  end

  @_help[:cmd] = "execute command within the target build dir, e.g. rake localGCC_cmd['pwd']" unless @_help.has_key?(:cmd)
  task toDo[:cmd] ,:cmd do |t, args|
    warn "No command given!!" && exit(1) if args.cmd.nil?
    execute(args.cmd,builder)
  end

  @_help[:mods] = "get the auto loaded modules on the target machine"
  task toDo[:mods] do |t|
    execute("module list", builder)
  end

  @_help[:log] = "log everything to <host|localTask>.log"
  file baseTaskName+".log" do |t|
    sh "rake #{baseTaskName} > #{t.name} 2>&1"
    sh "finished".colorize(color: :green)
  end

  # this is the mail task to be executed: sync, conf, make, check
  desc builder.docstring
  task baseTaskName.to_sym  => [syncSource ? toDo[:sync] : nil,
                                toDo[:conf],
                                toDo[:make],
                                toDo[:check]].compact.map(&:to_sym)
end
# }}}
def getUsername(builderConfig, hostConfig)
  username = nil

  username = @user if [builderConfig['hostname'],hostConfig['hostname']].include?('localhost')
  return username unless username.nil?

  username =  builderConfig['username'] if builderConfig.has_key?('username')
  return username unless username.nil?

  username = hostConfig['username'] if hostConfig.has_key?('username')
  return username unless username.nil?

  username = @userConfig["remoteUser"] if @userConfig.has_key?('remoteUser')
  return username unless username.nil?

  warn "Could not find username!!"
  exit(1)
end
# constuct builders out of user configuration {{{ ==============================
Builder = Struct.new(:host,:hostname,:username,:compiler,:targetDir,:configureCall,:isLocal?,:docstring,:envConfigFiles)
# 1) construct builders from host configuration
#    this is what config/default should be able to handle
#    CC must be given, otherwise it's just a host description
@userConfig["hosts"].each {|host,config|
  config['CC'].each {|cc|
    builder = Builder.new(host,
                          config["hostname"],
                          ('localhost' == config['hostname']) \
                              ? @user \
                              : ( config.has_key?('username') \
                                 ? config['username'] \
                                 : @userConfig["remoteUser"]),
                          cc,
                          [config["dir"],cc,getBranchName].join(File::SEPARATOR),
                          "./config/default CC=#{cc}",
                          config["hostname"] == 'localhost',
                         "builder on #{config['hostname']}, CC=#{cc}",
                          config.has_key?("envConfigFiles") ? config["envConfigFiles"] : [])

    builder2task(builder)
  } if config.has_key?('CC')
}
# 2) construct builders from manual configuration
@userConfig["builders"].each {|builderName,builderConfig|
  hostConfig = @userConfig['hosts'][builderConfig['hostname']]

  warn "Hostconfig not found!!!!" and exit(1) if hostConfig.nil?

  username, hostname = getUsername(builderConfig,hostConfig), hostConfig['hostname']

  if username.nil? or hostname.nil? then
    puts [username, hostname].join(' - ').colorize(color: :red)
    warn "Missing connection info!"
    exit(1)
  else
    puts [username, hostname].join(' - ').colorize(color: :green) if false
  end

  builder = Builder.new(builderName,
                        hostname,
                        username,
                        '', # CC can be empty here, because it should be set by the given configureCall
                        [hostConfig['dir'],builderName,getBranchName].join(File::SEPARATOR),
                        builderConfig['configureCall'],
                        [builderConfig['hostname'],hostConfig['hostname']].include?('localhost'),
                        builderConfig.has_key?('docstring') \
                          ? builderConfig['docstring'] \
                          : "builder on #{builderConfig['hostname']}: #{builderConfig['configureCall']}",
                        builderConfig.has_key?("envConfigFiles") ? builderConfig["envConfigFiles"] : [])

  builder2task(builder,true, builderConfig['sync'])

} if @userConfig.has_key?('builders')
# }}}
#

desc "execute listed tasks in parallel, each of them in a separate xterm"
task :par do |t|
  # remove all tasks from the stack
  Rake.application.top_level_tasks.clear

  # create a task list from the command line
  ARGV.shift
  taskList = ARGV
  dbg(taskList)

  # execute tasks in parallel
  Parallel.map(taskList) {|t| sh "xterm -hold -e 'rake #{t}' " }
end

desc "show help on all hidden tasks"
task :help do
  @_help.each {|t,help|
    sep = :log == t ? '.' : '_'
    puts "rake <host|localTask>#{sep}#{t}".ljust(35,' ') + "# #{help}"
  }
end

task :default do |t|
  sh "rake -sT"
end

desc "generate tags data base for vim amd emacs"
task :tags do |t|
  srcFiles = Dir.glob("src/**/*.{h,c}") + Dir.glob("libcdi/**/*.{c,h,cpp,hpp,f90,f}")
  Parallel.map(["","-e"]) {|ctagsOutputMode|
    sh "ctags #{ctagsOutputMode} #{srcFiles.join(' ')}"
  }
end
# check connections {{{
task :checkConnections do |t|
  pp Parallel.map(@userConfig["hosts"]) {|host, config|
    hostname = config['hostname']
    username = 'localhost' == config["hostname"] \
      ? @user \
      : ( config.has_key?('username') \
         ? config['username'] \
         : @userConfig["remoteUser"] )

    if 'localhost' == config["hostname"] then
      config["hostname"]
    else
      Net::SSH.start(hostname,username) do |ssh|
        hn = ssh.exec!("hostname -f").strip
        us = ssh.exec!("echo $USER").strip
        [hn,us]
      end
    end
  }
end
# }}}
