import gym
import numpy as np
from gym import error, spaces, utils
from gym.utils import seeding
import oscillator_cpp
from gym import Env
from gym.spaces import Discrete, MultiDiscrete, MultiBinary, Box


class oscillatorEnv(gym.Env):
  metadata = {'render.modes': ['human']}
  def __init__(self, len_state=150, ep_length=5000, nosc=1000, epsilon=0.2, frrms=0.02,ndim=3,random=False,
               sigmoid=False):

    """ 
    Init function:
    sigmoid: Function that we observe instead of original one: Bool
    len_state: shape of state that agent observes [250,1]: integer

    BVDP params
    nosc: number of oscillators: integer
    epsilon: coupling parameter: float
    frrms: width of the distribution of natural frequencies: float
    """
    
    super(oscillatorEnv, self).__init__()
    #print(1)
    #Call init function and save params

    self.y = oscillator_cpp.init(nosc,epsilon,frrms,ndim)
    #print(2)
    self.nosc = nosc
    self.epsilon = epsilon
    self.frrms = frrms
    self.ndim = ndim
    self.ep_length = ep_length

    #Dimensionality of our observation space
    self.dim = 1
    self.action_space = Box(low=-0.5, high=0.5, shape=(1,), dtype=np.float64)
    self.observation_space = Box(low=-2, high=2, shape=(len_state,), dtype=np.float64)

    #Meanfield for all neurons
    self.x_val = oscillator_cpp.Calc_mfx(self.y)
    self.y_val = oscillator_cpp.Calc_mfy(self.y)

    #Episode Done?
    self.done = False
    self.current_step = 0 

    #Our current state, with length(1,len_state)
    self.y_state = []
    self.x_state = []

    self.len_state = len_state
    
    #Reset environment
    self.reset()
    
  def step(self, action):
      """
      Function that called at each step.

      action: signal to make perturbation: [[float]]
      returns: arrayed_version:np.array(1,len_state), 
      Reward: Our reward function :float, 
      done: Does it end? :Bool, 
      additional_information: Nothing to show :( :{} 
      """
      #Vectorized form
      self.y_one = oscillator_cpp.get_one(self.y,1)
      self.y_900 = oscillator_cpp.get_one(self.y,25)
      val = float(action[0])
      if val != 0:
        self.y = oscillator_cpp.Pertrubation(self.y, val)
      self.y = oscillator_cpp.Make_step(self.y)
      
          
      #Calculate MeanField
      self.x_val = oscillator_cpp.Calc_mfx(self.y)
      self.y_val = oscillator_cpp.Calc_mfy(self.y)
      
      #if sigmoid:
          #self.x_val = sigmoid(self.x_val)
          #self.y_val = sigmoid(self.y_val)
      
      #Save our state
      self.y_state.append(self.y_val)
      self.x_state.append(self.x_val)
      
      #Check length of our state
      if len(self.y_state) > self.len_state:
          self.y_state = self.y_state[1:]
          self.x_state = self.x_state[1:]

      self.current_step += 1

      self.done = self.current_step >= self.ep_length

      #Make vectorized form
      arrayed_version = np.array(self.y_state)
      
      #if sigmoid:
          #arrayed_version = sigmoid(arrayed_version)
          
      return arrayed_version, self.Reward(self.x_val,np.array(self.x_state),self.y_val,val), self.done, {} 

  
  def reset(self):
    """
    Reset environment

    Returns:arrayed_version:np.array(1,len_state)
2
    """
    self.current_step = 0 
    self.y_state = []
    self.x_state = []
    self.y = oscillator_cpp.init(self.nosc,self.epsilon,self.frrms,self.ndim)
    for i in range(self.len_state):
        oscillator_cpp.Make_step(self.y)
        
        self.x_val = oscillator_cpp.Calc_mfx(self.y)
        self.y_val = oscillator_cpp.Calc_mfy(self.y)

        self.y_state.append(self.y_val)
        self.x_state.append(self.x_val)
        
        #Check length of our state
        if len(self.y_state) > self.len_state:
            self.y_state = self.y_state[1:]
            self.x_state = self.x_state[1:]

    arrayed_version = np.array(self.y_state)    
    #if sigmoid:
        #arrayed_version = sigmoid(arrayed_version)
    return arrayed_version
    
  def render(self, mode='human', close=False):
    """
    Pass...
    """
    pass

  def Reward(self, x,x_state,y,action_val):
    """
    Super duper reward function, i am joking, just sum of absolute values which we supress...
    returns: float
    """
    #Bursting
    #print(x)
    #print(-(x-np.mean(x_state))**2 - 2*np.abs(action_val))
    return -((x+0.4)**2)*2 - 0.5*np.abs(action_val) 
def sigmoid(x):
    x_0 = 1.2
    return 1. / (1. + np.exp(-x/x_0))

