"""
Created on Mar 7, 2017

@author: Mirela Andronescu (mandronescu@bccrc.ca)
"""
from kronos.utils import ComponentAbstract
import os


class Component(ComponentAbstract):
    '''
    eval_epiclomal
    '''

    def __init__(self, component_name='eval_epiclomal', component_parent_dir=None, seed_dir=None):
        self.version = "1.0.0"

        ## initialize ComponentAbstract
        super(Component, self).__init__(component_name, component_parent_dir, seed_dir)

    def make_cmd(self, chunk=None):

        cmd = self.requirements['Rscript']
        seed = os.path.join(self.seed_dir, "eval_epiclomal.R")
        hdist_software = self.requirements['hdist_software']
        visualization_software = self.requirements['visualization_software']

        # cmd_args = [
            # seed,
            # self.args.input_dir,
            # self.args.output_dir,
            # self.args.model, 
            # vmeasure_software,
            # hdist_software
        # ]
        
        cmd_args = [seed]

        for k, v in vars(self.args).iteritems():
            if v is None or v is False:
                continue
            cmd_args.append("--" + k + " " + str(v))
            
        cmd_args.append("--hdist_software " + hdist_software)    
        cmd_args.append("--visualization_software " + visualization_software)         
        return cmd, cmd_args

## to run as stand alone
def _main():
    eval_epiclomal = Component()
    eval_epiclomal.args = component_ui.args
    eval_epiclomal.run()


if __name__ == '__main__':
    import component_ui

    _main()
