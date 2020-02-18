import gym
import numpy as np
from gym import error, spaces, utils
from gym.utils import seeding
import oscillator_cpp
from gym import Env
from gym.spaces import Discrete, MultiDiscrete, MultiBinary, Box
from stable_baselines.common.policies import MlpPolicy,MlpLnLstmPolicy,FeedForwardPolicy
from stable_baselines.common.vec_env import DummyVecEnv,SubprocVecEnv,VecNormalize, VecEnv
from stable_baselines import PPO2
from stable_baselines.common.vec_env import VecEnv

class oscillatorEnv(gym.Env):
  metadata = {'render.modes': ['human']}
  def __init__(self, len_state=250, ep_length=10000, nosc=1000, epsilon=0.03, frrms=0.1,random=False,
                model = None,
               initial_steps=55000,
                model_steps=55000,
                skip_rate=0,
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
    #Call init function and save params
    if random:
        epsilon = np.random.uniform(0.03,0.5)
    self.y = oscillator_cpp.init(nosc,epsilon,frrms)
    self.nosc = nosc
    self.epsilon = epsilon
    self.frrms = frrms
    
    self.ep_length = ep_length
    self.model_steps = model_steps
    #Dimensionality of our observation space
    self.dim = 1
    self.action_space = Box(low=-1, high=1, shape=(1,), dtype=np.float32)
    self.observation_space = Box(low=-1.5, high=1.5, shape=(len_state,), dtype=np.float32)

    #Meanfield for all neurons
    self.x_val = oscillator_cpp.Calc_mfx(self.y)
    self.y_val = oscillator_cpp.Calc_mfy(self.y)

    #Episode Done?
    self.done = False
    self.current_step = 0 

    #Our current state, with length(1,len_state)
    self.y_state = []
    self.x_state = []
    #Our actions 
    self.actions = []
    

    self.skip_rate = skip_rate
    self.skip = 1


    self.len_state = len_state
    self.initial_steps = initial_steps
    self.model = model
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
           
      if self.skip_rate >= 1:
          val = 0
          if self.skip > self.skip_rate:  
            val = float(action[0])
            self.skip = 1
          else:
            self.skip += 1
          self.val = val
      else:
        val = float(action[0])
      self.actions.append(val)
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
        
      self.x_states.append(self.x_val)
    
      #Check length of our state
      if len(self.y_state) > 250:
          self.y_state = self.y_state[1:]
          self.x_state = self.x_state[1:]

      self.current_step += 1

      self.done = self.current_step >= self.ep_length

      #Make vectorized form
      arrayed_version = np.array(self.y_state)

      return arrayed_version, self.Reward(self.x_val,self.x_state,val), self.done, {} 

  def reset(self):
    """
    Reset environment

    Returns:arrayed_version:np.array(1,len_state)

    """
    
    self.current_step = 0 
    self.y_state = []
    self.x_state = []
    self.x_states=[]
    self.actions = []
    self.y = oscillator_cpp.init(self.nosc,self.epsilon,self.frrms)
    self.skip = 1
    for i in range(self.initial_steps):
        oscillator_cpp.Make_step(self.y)
        
        self.x_val = oscillator_cpp.Calc_mfx(self.y)
        self.y_val = oscillator_cpp.Calc_mfy(self.y)
        self.actions.append(0)
        self.y_state.append(self.y_val)
        self.x_state.append(self.x_val)
        
        self.x_states.append(self.x_val)
        
        #Check length of our state
        if len(self.y_state) > 250:
            self.y_state = self.y_state[1:]
            self.x_state = self.x_state[1:]
        
    
    if self.model != None:
        for i in range(self.model_steps):
            action, _states = self.model.predict( np.array(self.y_state))
            self.actions.append(action[0])
            val = float(action[0])
            self.y = oscillator_cpp.Pertrubation(self.y, val)
            self.y = oscillator_cpp.Make_step(self.y)


            #Calculate MeanField
            self.x_val = oscillator_cpp.Calc_mfx(self.y)
            self.y_val = oscillator_cpp.Calc_mfy(self.y)

            #if sigmoid:
              #self.x_val = sigmoid(self.x_val)
              #self.y_val = sigmoid(self.y_val)

           
            self.y_state.append(self.y_val)
            self.x_state.append(self.x_val)
            
            self.x_states.append(self.x_val)

            if len(self.y_state) > 250:
                self.y_state = self.y_state[1:]
                self.x_state = self.x_state[1:]
            

    arrayed_version = np.array(self.y_state)
        
    return arrayed_version
    
  def render(self, mode='human', close=False):
    """
    Pass...
    """
    pass

  def Reward(self, x,x_state,action_value):
    """
    Super duper reward function, i am joking, just sum of absolute values which we supress...
    returns: float
    """
    #print(baseline)
    #print(-(x-np.mean(x_state))**2 - 5*np.abs(action_value))
    return -(x-np.mean(x_state))**2 - 2*np.abs(action_value)

