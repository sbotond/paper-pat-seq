import      os
import      sys
import      time
import      copy
import      subprocess                      as      sp
from        Bio                             import  SeqIO
import      matplotlib                                                                                                                                         
matplotlib.use('PDF')
from        matplotlib.backends.backend_pdf import  PdfPages
from        matplotlib                      import  pyplot          as  plt
from        collections                     import  defaultdict
import      numpy                           as      np
import      cPickle
import      scipy.stats                     as      st

class LinReg:
    """ Linear regression """
    def __init__(self, x, y):
        self.x      = x
        self.y      = y

        self.slope      = None
        self.intercept  = None
        self.r          = None
        self.p          = None
        self.se         = None

        self.reg()

    def reg(self):
        """ Perform linear regression """
        (self.slope, self.intercept, self.r, self.p, self.se) = st.linregress(self.x, self.y)

    def trendline(self,minx,maxx,col='blue'):
        """ Add trendline """
        x   = np.linspace(minx, maxx)
        y   = self.intercept + self.slope * x
        plt.plot(x,y,'-',color=col)

class Report:
    """ Class for plotting reports """
    def __init__(self, pdf):
        self.pdf    = pdf
        self.pages  = PdfPages(pdf)

    def plot_tails(self, tr, nv, gv, g_runs):
        title=tr
        xlab="Position"
        ylab="Fragment coverage"
        la="NVTR"
        lb="G-tailed"

        fig = plt.figure()
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212)
        fig.set_tight_layout(True)

        x   = np.arange(len(nv))
        ax1.plot(x, nv, 'b.', label=la)
        ax1.plot(x, gv, 'g^', label=lb)
        ax1.set_xlim((0,len(x)))
        mc  = max(np.max(nv),np.max(gv))
        ax1.set_ylim((0.98,mc+1))
        ax1.set_xlabel(xlab)
        ax1.set_ylabel(ylab)
        ax1.legend()
        ax1.set_title(title)
        
        # ax1.clf()

        if tr in g_runs:
            title="Tail run distribution"
            xlab="Length"
            ylab="Count"

            ax2.bar(g_runs[tr].keys(), g_runs[tr].values(), width=0.1)
            ax2.set_xlabel(xlab)
            ax2.set_ylabel(ylab)
            ax2.set_title(title)

            # ax2.clf()

        plt.close(fig)
        self.pages.savefig(fig)

    def plot_hash(self, h, title="", xlab="", ylab=""):
        """ Visualise hash as a bar plot """
        fig = plt.figure()
        plt.bar(h.keys(), h.values(), width=0.1)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def finalise_plot(self, fig):
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def plot_hist_core(self, a, title="", xlab="", ylab="",bins=50,col="blue",alpha=0.5,mark_mean=True):
        """ Plot histogram """
        fig = plt.figure()
        mean  = np.mean(a)
        plt.hist(a,bins=bins,color=col,alpha=alpha,label="Mean: %g" % mean)
        plt.axvline(mean,color=col) 
        plt.legend()
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)

        return fig

    def plot_hist(self, *args, **kwargs):
        fig = self.plot_hist_core(*args, **kwargs)
        self.finalise_plot(fig)

    def plot_hist_with_readlen(self, *args, **kwargs):
        run_len = kwargs.pop("run_len")
        fig = self.plot_hist_core(*args, **kwargs)
        plt.axvline(run_len, color="red")
        self.finalise_plot(fig)

    def plot_double_hists(self,a,b,title="",lab1="1",lab2="2",xlab="",ylab="",bins=100,acol="blue",bcol="red",alpha=0.5,stacked=False,center_zero=False,normed=False,mark_means=False,force_int=True):
        """ Overlayed plot of two histograms """
        a=np.array(a)
        b=np.array(b)
        fig = plt.figure()

        if center_zero:
            # Set 0.0 as a boundary:
            u = max( np.max(a), np.max(b))
            l = min( np.min(a), np.min(b))
            bp = np.arange(0,u,2.0/bins)
            bn = np.arange(0,l,-2.0/bins)[1:][::-1]
            bins  = np.hstack((bn,bp))

        if force_int:
            u = int(max( np.max(a), np.max(b)))
            l = int(min( np.min(a), np.min(b)))
            bins=range(l,u)

        plt.hold(True)
        if len(a) > 0:
            amed   = np.mean(a)
            plt.hist(a, bins=bins,color=acol, stacked=stacked,alpha=alpha, label="%s (mean: %g)" % (lab1, amed), align='mid',normed=normed)
            if mark_means:
               plt.axvline(amed,color=acol) 
        if len(b) > 0:
            bmed   = np.mean(b)
            plt.hist(b, bins=bins,color=bcol, stacked=stacked,alpha=alpha, label="%s (mean: %g)" % (lab2, bmed), align='mid',normed=normed)
            if mark_means:
               plt.axvline(bmed,color=bcol) 
        plt.hold(False)

        plt.legend()
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def plot_contour(self, z, title="", xlab="", ylab=""):
        """ Visualise matrix as a filled contour plot """
        fig = plt.figure()
        p   = plt.contourf(z)
        plt.colorbar(p, orientation='vertical')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def plot_array(self, y, title="", xlab="", ylab=""):
        """ Visualise  array as a bar plot """
        fig = plt.figure()
        plt.bar(np.arange(len(y)),y,width=0.1)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)


    def plot_arrays(self,x, y, title="", xlab="", ylab="",style='.',col='red',alpha=0.5,trendline=True):
        """ Visualise  array as a bar plot """
        fig = plt.figure()
        plt.hold(True)
        plt.plot(x,y,style,color=col,alpha=alpha)
        lr  = LinReg(x,y)
        lr.trendline(minx=np.min(x), maxx=np.max(x), col=col)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.suptitle(title)
        plt.title("r = %g p-value = %g se = %g" % (lr.r, lr.p, lr.se) )
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)

    def plot_liks(self,x,la,lb,lj,na="",nb="",nj="",title="", xlab="", ylab=""):
        """ Visualise likelihood profiles """
        fig = plt.figure()
        plt.hold(True)
        cut = min(max(la), max(lb), max(lj))
        plt.ylim((cut * 5,0))
        plt.plot(x,la,'-',label=na)
        plt.plot(x,lb,'-',label=nb)
        plt.plot(x,lj,'-',label=nj)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.legend()
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)


    def plot_arrays_double(self,a, b, title="", xlab="", ylab="", la="", lb=""):
        """ Visualise two array as overlapping bar plots """
        fig = plt.figure()
        x   = np.arange(len(a))
        plt.plot(x, a, 'b.', label=la)
        plt.plot(x, b, 'g^', label=lb)
        plt.xlim((0,len(x)))
        mc  = max(np.max(a),np.max(b))
        plt.ylim((0.98,mc+1))
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.legend()
        plt.title(title)
        self.pages.savefig(fig)
        plt.clf()
        plt.close(fig)


    def close(self):
        self.pages.close()

class Log:
    """ Logging utility class """
    def __init__(self, fname=None, level=0):
        self.level = level
        if fname == None:
            self.fname  = "<sys.stdout>"     
            self.file   = sys.stdout
        else:
            self.file   = open(fname, "w")
            self.fname  = fname

    def close(self):
        self.file.flush()
        self.file.close()

    def log(self, message):
        if self.level < 0:
            return
        self.file.write("[%s] %s\n" % (time.strftime("%y-%m-%d %H:%M:%s"), message) )

    def fatal(self, message):
        self.file.write("[%s] %s\n" % (time.strftime("%y-%m-%d %H:%M:%s"), message) )
        sys.exit(1)

class Rtemp:
    """ Utility class handling temporary storage """
    def __init__(self, base, log, autoclean=False):
        self.log        = log
        self.autoclean  = autoclean
        self.parent     = None
        self.children   = [ ]
        self.files      = [ ]
        if os.path.isdir(base) != True:
            log.fatal("The base must be a directory: %s" % base)
        self.base   = os.path.abspath(base)

    def exists(self, fname):
        """ Check wheteher a file exists """
        return os.path.exists(fname) 

    def _iterate_fname(self, fname):
        """ Iterate until we don't have a name clash """
        i           = 0
        orig_fn     = fname # basename
        while True:
            if (fname in self.files) or self.exists(fname):
                i       += 1
                fname   = orig_fn + ("_%03d" % i)
            else:
                break
        return fname

    def tempfile(self, name):
        """ Get a temporary file """
        fname       = self._iterate_fname(os.path.join(self.base, name))
        self.register(fname)            
        return fname

    def temp_fh(self, name):
        """ Get a temporary file handle """
        fname   = self.tempfile(name)
        f       = open(fname, "w")
        return f

    def clean(self):
        """ Remove registered temporary files """
        for child in self.children:
            child.clean()   # call cleanup on children.
        tmp = list(self.files)
        for f in tmp:
            self.remove(f)
        # Delete the directory if children is  a subdir:
        if self.parent != None:
            os.rmdir(self.base)

    def remove(self, fname):
        """ Remove a temporary file """
        if not (fname in self.files):
            self.log.fatal("The file %s is not mannaged by this object!" % fname)
        if os.path.exists(fname):
            os.remove(fname)
        self.files.remove(fname)

    def subdir(self, dname):
        """ Get a mannaged temporary subdirectory """
        clone           = copy.copy(self)        
        dname           = os.path.join(self.base, dname)
        clone_base      = self._iterate_fname(dname)
        os.mkdir(clone_base)
        clone.base      = clone_base
        clone.parent    = self
        clone.children  = []
        clone.files     = []
        self.children.append(clone)
        return clone

    def register(self, fname):
        """ Register temporary file """
        self.files.append(fname)
    
    def unregister(self, fname):
        """ Unregister temporary file """
        self.files.remove(fname)

    def __del__(self):
        if self.autoclean:
            self.clean()

class Cmd:
    """ Class for running commands """
    def __init__(self, log, prog, pre_args=[], post_args=[], opts={}, inp=None, outp=None, cwd=None, files=[ ], sink_err=False, autoclean=False, path=None):
        self.log        =   log         # log object.
        self.pre_args   =   pre_args    # before options arguments
        self.post_args  =   post_args   # after options arguments
        self.opts       =   opts        # option flags
        self.cwd        =   cwd         # execution directory
        self.files      =   files       # registered files
        self.sink_err   =   sink_err    # sink error messages
        self.autoclean  =   autoclean   # automatic cleanup
        self.path       =   path        # path to binary
        self.name       =   os.path.basename(prog)
        if self.path == None:
            self.prog       =   prog
        else:
            self.prog       =   os.path.join(path, prog)
        self.pre_args   =   pre_args    # before options arguments

        # Build command:
        self.build_cmd()

        # Sort out input:
        input_fh    = sp.PIPE   # get input from pipe
        if inp != None:
            if type(inp) == file:
                input_fh    = inp   # input is a file handler
            elif type(inp) == str:
                input_fh    = open(inp, "r") # input is a file name
            else:
                self.log.fatal("Invalid input source!")
        self.input_fh   = input_fh

        # Sort out output:
        output_fh = sp.PIPE # Output to pipe
        if outp != None:
            if type(outp) == file:
                output_fh = outp    # file handler output
            elif type(outp) == str:
                output_fh    = open(outp, "w")  # create file
            else:
                self.log.fatal("Invalid output!")
        self.output_fh  = output_fh

        # Open pipe:
        self.pipe   =   sp.Popen(args=self.cmd,bufsize=0,stderr=sp.PIPE,stdin=input_fh, stdout=output_fh, shell=False,cwd=self.cwd)

    def build_cmd(self):
        """ Build command list """
        cmd  = [ self.prog ] # executable
        cmd.extend(self.pre_args) # pre-arguments
        # Append flags:
        for (flag, val) in self.opts.iteritems():
            # Flag:
            cmd.append(flag)
            # Value if not none:
            if val != None:
                cmd.append(val)
        # Post-args:
        cmd.extend(self.post_args)
        self.cmd    = cmd

    def comm(self, in_data=''):
        """ Communicate with subprocess """
        (stdout_data, stderr_data)  = self.pipe.communicate(in_data)
        # Log stderr output:
        if stderr_data != None and self.sink_err != True:
            self.log_data(stderr_data)
        return stdout_data

    def log_data(self, data):
        """ Logging helper """
        log = self.log
        tmp = data.split('\n')  
        for line in tmp:
           log.log("<%s> %s" % (self.name, line)) 

    def register(self, fn):
        """ Register temporary file """
        self.files.append(fn)

    def unregister(self, fn):
        """ Unregister temporary file """
        self.files.remove(fn)

    def clean(self):
        """ Cleanup registered temporary files """
        for f in self.files:
            if os.path.exists(f):
                os.remove(f)
        self.files  = [ ]

    def __del__(self):
        if (self.pipe.stdin != None) and (not self.pipe.stdin.closed):
            # Close stdin and wait for final output.
            self.pipe.stdin.close()
            self.pipe.wait()
        if self.autoclean:
            self.clean()

    def close(self):
        # Close stdin, wait for final output.
        if self.pipe.stdin != None:
            self.pipe.stdin.close()
        self.pipe.wait()

class Pipe:
    """ Pipe helper class """
    def __init__(self):
        self.r, self.w  = os.pipe()

class Fasta:
    """ Fasta parsing class """
    def __init__(self, infile):
        self.infile     = infile
        self.in_fh      = open(infile, "r")
        self.iter       = SeqIO.parse(self.in_fh,'fasta')

    def __iter__(self):
        """ Return iterator """
        return iter(self.iter)

    def slurp(self):
        """ Slurp sequences """
        records = { }
        for s in iter(self):
            records[s.name] = str(s.seq)
        return records

def chext(s, new_ext):
        """ Change extension of a file name """
        tmp = s.split('.')
        tmp = tmp[0:-1]
        tmp.append(new_ext)
        return '.'.join(tmp)

def unpickle(f):
    """ Load data from pickle file """
    fh     = open(f, "r")
    pickle = cPickle.Unpickler( fh )
    tmp = pickle.load()
    return tmp

def pickle(d, fname):
    """ Pickle data to file """
    fh     = open(fname, "w")
    pickle = cPickle.Pickler( fh )
    pickle.dump(d)
    fh.flush()
    fh.close()

class TranscriptList:
    """ Transcript white list """
    def __init__(self, infile, log):
        self.log        = log
        self.have_list  = False
        if not infile is None:
            self.infile = infile
            fh          = open(infile, "r")
            self.trs    = fh.readlines()
            for i in xrange(len(self.trs)):
                self.trs[i] = self.trs[i].rstrip()
            self.have_list = True
            self.log.log("Number of target transcripts specified: %d" % self.size() )
        

    def in_list(self, name):
        if not self.have_list:
            return True
        else:
            return (name in self.trs) 

    def size(self):
        return len(self.trs)

