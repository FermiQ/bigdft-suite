"""Calculation datasets.

This module deals with the handling of series of calculations.
Classes and functions of this module are meant to simplify the approach to ensemble calculations
with the code, and to deal with parallel executions of multiple instances of the code.

"""

class Dataset(Runner):
    """A set of calculations.

    Such class contains the various instances of a set of calculations with the code.
    The different calculations are labelled by parameter values and information that might then be
    retrieved for inspection and plotting.
    """
    def __init__(self,label='BigDFT dataset',run_dir='runs',input={}):
        """
        Set the dataset ready for appending new runs

        Args:
           label (str): The label of the dataset. It will be needed to identify it in plotting 
           run_dir (str): path of the directory where the runs will be performed.
           input (dict): Inputfile to be used for the runs, can be overridden by the inputs of the run
        """
        from copy import deepcopy
        from futile.Utils import make_dict
        Runner.__init__(self,label=label,run_dir=run_dir,input=make_dict(input))
        self.runs=[]
        """List of the runs which have to be treated by the dataset. Set is organized by calculators in order to
           run separate instances of the calculators to be performed
        """
        self.calculators=[]
        """ 
        Calculators which will be used by the run method, useful to gather the inputs in the case of a multiple sending.
        """
        self.results={}
        """ Set of the results of each of the runs
        """
    def append_run(self,id,runner,**kwargs):
        """Add a run into the dataset.
        
        Append to the list of runs to be performed the corresponding runner and the arguments which are associated to it.
     
        Args:
          id (dict): the id of the run, useful to identify the run in the dataset. It has to be a dictionary as it may contain
             different keyword. For example a run might be classified as ``id = {'hgrid':0.35, 'crmult': 5}``.
          runner (Runner): the runner class to which the remaining keyword arguments will be passed at the input.
        """
        from copy import deepcopy
        #get the number of this run
        irun=len(self.runs)
        #create the input file for the run, combining run_dict and input
        inp_to_append=deepcopy(kwargs)
        inp_to_append.update(self.global_options['input'])
        inp_to_append['run_dir']=self.global_options['run_dir']
        #append it to the runs list
        self.runs.append(inp_to_append)
        #search if the calculator already exists
        found = False
        for calc in self.calculators:
            if calc['calc'] == runner:
                calc['runs'].append[irun]
                found=True
                break
        if not found:
            self.calculators.append({'calc': runner, 'runs':[irun]})

    def process_run(self):
        """
        Run the dataset, by performing explicit run of each of the item of the runs_list.
        """
        for c in self.calculators:
            calc=c['calc']
            #we must here differentiate between a taskgroup run and a separate run
            for r in c['runs']:
                inp=self.inputs[r]
                self.results[r]=calc.run(**inp)
        return {}
    
            

