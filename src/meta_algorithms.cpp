#include "utils.h"
#include "ucb.h"
#include "meta_algorithms.h"
#include <fstream>
#include <algorithm>    // std::sort

MetaAlgorithm::MetaAlgorithm(string name, int num_of_arms) : MABAlgorithm(name, num_of_arms) {}

MABAlgorithm* MetaAlgorithm::get_sub_alg() {
  MetaAlgorithm* b = dynamic_cast<MetaAlgorithm*>(this->sub_alg);
  if (b == NULL) { // if not meta
    return this->sub_alg;
  } else { // if meta
    return b->get_sub_alg();
  }
}

Round_Algorithm::Round_Algorithm(string name, int num_of_arms, string sub_alg_line, boost::mt19937& rng, MAB* mab) : MetaAlgorithm(name, num_of_arms) {
  this->sub_alg = get_algorithm(sub_alg_line, &rng, mab);
  this->reset(-1);
}

void Round_Algorithm::reset(int action) {
	this->MABAlgorithm::reset(action);
  if (action == -1) {
    this->arms_pulled.assign(this->num_of_arms, false);
  } else {
    this->arms_pulled[action] = false;
  }
  this->sub_alg->reset(action);
}

int Round_Algorithm::choose_action() {
	int arm_to_pull = -1;
  for (int i = 0; i < this->arms_pulled.size(); i++) {
    if (!this->arms_pulled[i]) {
      arm_to_pull = i;
      break;
    }
  }
	if (arm_to_pull == -1) {
		// Delegate to sub_alg
		arm_to_pull = this->sub_alg->choose_action();
	}
	return arm_to_pull;
}

void Round_Algorithm::receive_reward(double reward, int pulled_arm) {
  this->MABAlgorithm::receive_reward(reward, pulled_arm);
  this->arms_pulled[pulled_arm] = true;
	this->sub_alg->receive_reward(reward, pulled_arm);
}




Algorithm_With_Uniform_Exploration::Algorithm_With_Uniform_Exploration(string name, int num_of_arms, string sub_alg_line, double alpha, boost::mt19937& rng, MAB* mab) : MetaAlgorithm(name, num_of_arms) {
	this->sub_alg = get_algorithm(sub_alg_line, &rng, mab);
	this->alpha = alpha;
}

void Algorithm_With_Uniform_Exploration::reset(int action) {
	this->MABAlgorithm::reset(action);
  this->sub_alg->reset(action);
}

int Algorithm_With_Uniform_Exploration::choose_action() {
	int arm_to_pull = -1;
	float r = random_unit();
	if (r < this->alpha) { // explore
		arm_to_pull = rand() % this->num_of_arms;
	} else {
		arm_to_pull = this->sub_alg->choose_action();
	}
	return arm_to_pull;
}

void Algorithm_With_Uniform_Exploration::receive_reward(double reward, int pulled_arm) {
  this->MABAlgorithm::receive_reward(reward, pulled_arm);
	this->sub_alg->receive_reward(reward, pulled_arm);
}



CD_Algorithm::CD_Algorithm(string name, int num_of_arms, int M, string cdt_line, string sub_alg_line, bool use_history, int max_history, bool smart_resets, boost::mt19937& rng, MAB* mab) : MetaAlgorithm(name, num_of_arms) {
  this->M = M;
  this->max_history = max_history;
  this->sub_alg = get_algorithm(sub_alg_line, &rng, mab);
  for (int i = 0; i < num_of_arms; i++) {
    this->cdts_up.push_back(get_cdt(cdt_line + " 1"));
    this->cdts_down.push_back(get_cdt(cdt_line + " 0"));
    if (smart_resets) {
      this->sr_cdts_up.push_back(get_cdt(cdt_line + " 1"));
      this->sr_cdts_down.push_back(get_cdt(cdt_line + " 0"));
    }
  }
  this->use_history = use_history;
  this->smart_resets = smart_resets;
  if (smart_resets) {
    this->sr_change_estimates.assign(num_of_arms, 0);
  }
  this->last_resets_up.assign(num_of_arms, false);
  //this->num_alarms_up = 0;
  //this->num_alarms_down = 0;
  this->reset(-1);
}

void CD_Algorithm::reset(int action) {
  this->MABAlgorithm::reset(action);
  this->sub_alg->reset(action);
  if (action == -1) {
    for (int i = 0; i < this->num_of_arms; i++) {
      this->cdts_up[i]->reset(0);
      this->cdts_down[i]->reset(0);
      if (this->smart_resets) {
        this->sr_cdts_up[i]->reset(0);
        this->sr_cdts_down[i]->reset(0);
      }
    }
    this->till_M.clear();
    this->till_M.assign(this->num_of_arms, this->M);
    this->collected_rewards.clear();
    this->collected_rewards.resize(this->num_of_arms);
    if (this->smart_resets) {
      this->sr_change_estimates.assign(num_of_arms, 0);
      this->sr_collected_rewards.clear();
      this->sr_collected_rewards.resize(this->num_of_arms);
    }
    this->last_resets_up.assign(num_of_arms, false);
    this->timestep = 0;
    this->to_be_reset.clear();
  } else {
    this->cdts_up[action]->reset(0);
    this->cdts_down[action]->reset(0);
    this->till_M[action] = this->M;
    this->collected_rewards[action].clear();
    if (this->smart_resets) {
      this->sr_cdts_up[action]->reset(0);
      this->sr_cdts_down[action]->reset(0);
      this->sr_change_estimates[action] = 0;
      this->sr_collected_rewards[action].clear();
    }
    this->last_resets_up[action] = false;
    this->to_be_reset.erase(action);
  }
  //cout << "Alarms up: " << this->num_alarms_up / 100.0 << endl;
  //cout << "Alarms down: " << this->num_alarms_down / 100.0 << endl;
}

int CD_Algorithm::choose_action() {
  int action = -1;
  for (int i = 0; i < this->num_of_arms; i++) {
    if (this->till_M[i] > 0) {
      action = i;
      this->till_M[i]--;
      break;
    }
  }
  if (action == -1) {
    action = this->sub_alg->choose_action();
  }
  return action;
}

void CD_Algorithm::reset_arm(int pulled_arm) {
  this->sub_alg->reset(pulled_arm);
  this->MABAlgorithm::reset(pulled_arm);
  if (this->use_history) {
    int change_estimate = 0;
    bool alarm_up = this->last_resets_up[pulled_arm];
    //if (this->smart_resets) change_estimate = this->sr_change_estimates[pulled_arm];
    if (alarm_up) change_estimate = this->cdts_up[pulled_arm]->change_estimate;
    else change_estimate = this->cdts_down[pulled_arm]->change_estimate;

    /*if (this->smart_resets) {
      int history_amount = this->sr_collected_rewards[pulled_arm].size() - change_estimate;
      int t_start = change_estimate;
      if (history_amount > this->max_history) {
        t_start = this->sr_collected_rewards[pulled_arm].size() - this->max_history;
      }
      for (int i = t_start; i < this->sr_collected_rewards[pulled_arm].size(); i++) {
        this->sub_alg->receive_reward(this->sr_collected_rewards[pulled_arm][i], pulled_arm);
      }
    } else {*/
    int history_amount = this->collected_rewards[pulled_arm].size() - change_estimate;
    int t_start = change_estimate;
    if (history_amount > this->max_history) {
      t_start = this->collected_rewards[pulled_arm].size() - this->max_history;
    }
    for (int i = t_start; i < this->collected_rewards[pulled_arm].size(); i++) {
      this->sub_alg->receive_reward(this->collected_rewards[pulled_arm][i], pulled_arm);
    }
    //}
  }
  this->collected_rewards[pulled_arm].clear();
  //if (this->smart_resets) this->sr_collected_rewards[pulled_arm].clear();
  this->cdts_up[pulled_arm]->reset(0);
  this->cdts_down[pulled_arm]->reset(0);
  this->till_M[pulled_arm] = this->M;
}

void CD_Algorithm::receive_reward(double reward, int pulled_arm) {
  this->MABAlgorithm::receive_reward(reward, pulled_arm);
  this->sub_alg->receive_reward(reward, pulled_arm);
  this->collected_rewards[pulled_arm].push_back(reward);
  bool print_bounds = false;
  // // CODE FOR WRITING DOWN THE BOUNDS FOR THE MEANS.
  if(print_bounds){
    std::fstream outB("temp/exp12Bounds_"+this->name+".txt",std::ios::app);
    int tot_pulls = accumulate(this->num_of_pulls.begin(), this->num_of_pulls.end(), 0.);
    UCB1* u = dynamic_cast<UCB1*>(this->get_sub_alg());
    for(int a = 0; a != this->num_of_arms; ++a){
      double Q = u->means[a];
      double B = sqrt((2*log(tot_pulls + 1)) / max(this->num_of_pulls[a],1));
      outB << "arm " << a << " Q " << Q << " B " << B;
      outB << endl;
    }
  }
  

  bool alarm_up = this->cdts_up[pulled_arm]->run(reward);
  bool alarm_down = this->cdts_down[pulled_arm]->run(reward);
  //if (alarm_up) this->num_alarms_up++;
  //if (alarm_down) this->num_alarms_down++;
  bool alarm = alarm_up || alarm_down;
  if (alarm_up) this->last_resets_up[pulled_arm] = true;
  else if (alarm_down) this->last_resets_up[pulled_arm] = false;

  /*
  bool sr_alarm_up = false;
  bool sr_alarm_down = false;
  if (this->smart_resets) {
    this->sr_collected_rewards[pulled_arm].push_back(reward);
    sr_alarm_up = this->sr_cdts_up[pulled_arm]->run(reward);
    sr_alarm_down = this->sr_cdts_down[pulled_arm]->run(reward);
    if (sr_alarm_up) this->sr_change_estimates[pulled_arm] = this->sr_cdts_up[pulled_arm]->change_estimate;
    else if (sr_alarm_down) this->sr_change_estimates[pulled_arm] = this->sr_cdts_down[pulled_arm]->change_estimate;
  }*/

  if (alarm) {
    if (!this->smart_resets) {
      this->reset_arm(pulled_arm);
    } else {
      // Directly check UCB1...
      UCB1* u = dynamic_cast<UCB1*>(this->get_sub_alg());
      int best_action = -1;
      if (u != NULL) {
        best_action = distance(u->means.begin(), max_element(u->means.begin(), u->means.end()));
      }

      this->to_be_reset.insert(pulled_arm);

      if (pulled_arm == best_action && alarm_down) {
        for (auto a : this->to_be_reset) {
          this->reset_arm(a);
        }
        this->to_be_reset.clear();
      } else if (pulled_arm != best_action && alarm_up) {
        vector<int> arms_to_reset;
        arms_to_reset.push_back(pulled_arm);
        if (this->to_be_reset.find(best_action) != this->to_be_reset.end()) {
          arms_to_reset.push_back(best_action);
        }
        for (auto a : arms_to_reset) {
          this->reset_arm(a);
          this->to_be_reset.erase(a);
        }
      } else {
        this->cdts_up[pulled_arm]->reset(1);
        this->cdts_down[pulled_arm]->reset(1);
      }
    }
  }

  // sr cdts
  /*
  if (this->smart_resets && (sr_alarm_up || sr_alarm_down)) {
    this->sr_cdts_up[pulled_arm]->reset(0);
    this->sr_cdts_down[pulled_arm]->reset(0);
    this->sr_collected_rewards[pulled_arm].clear();
  }*/

  this->timestep++;
}

ADAPT_EVE::ADAPT_EVE(string name, int num_of_arms, int meta_duration, string cdt_line, string sub_alg_line, boost::mt19937& rng, MAB* mab) : MetaAlgorithm(name, num_of_arms) {

  this->meta_duration = meta_duration;
  this->rng = &rng;

  for (int i = 0; i < this->num_of_arms; i++) {
    this->cdts.push_back(get_cdt(cdt_line));
  }

  this->sub_alg = get_algorithm(sub_alg_line, this->rng, mab);
  this->other_sub_alg = get_algorithm(sub_alg_line, this->rng, mab);

  vector<Distribution*> fake_distributions;
  fake_distributions.push_back(new FixedDistribution("", 0, rng));
  fake_distributions.push_back(new FixedDistribution("", 0, rng));
  MAB* fake_mab = new MAB(fake_distributions, mab->timesteps);
  this->meta_alg = get_algorithm(sub_alg_line, this->rng, fake_mab);

  this->reset(-1);
}

void ADAPT_EVE::reset(int action) {
	this->MABAlgorithm::reset(action);
  if (action == -1) {
    for (int i = 0; i < this->num_of_arms; i++) {
      this->cdts[i]->reset(0);
    }
  } else {
    this->cdts[action]->reset(0);
  }

  this->sub_alg->reset(action);
  this->other_sub_alg->reset(action);
  this->meta_alg->reset(action);

  this->is_meta = false;
  this->t_meta = 0;

  this->timestep = 0;
}

int ADAPT_EVE::choose_action() {
  int action = -1;
  if (!this->is_meta) {
    action = this->sub_alg->choose_action();
  } else {
    this->saved_meta_action = this->meta_alg->choose_action();
    if (this->saved_meta_action == 0) {
      action = this->sub_alg->choose_action();
    } else {
      action = this->other_sub_alg->choose_action();
    }
  }
  return action;
}

void ADAPT_EVE::receive_reward(double reward, int pulled_arm) {
  this->MABAlgorithm::receive_reward(reward, pulled_arm);
  this->timestep++;

  if (!this->is_meta) {
    this->sub_alg->receive_reward(reward, pulled_arm);

    bool alarm = this->cdts[pulled_arm]->run(reward);
    if (alarm) {
      this->is_meta = true;
      this->t_meta = 0;

      this->other_sub_alg->reset(-1);
      this->meta_alg->reset(-1);
    }
  } else {
    t_meta++;

    this->meta_alg->receive_reward(reward, this->saved_meta_action);
    if (this->saved_meta_action == 0) {
      this->sub_alg->receive_reward(reward, pulled_arm);
    } else {
      this->other_sub_alg->receive_reward(reward, pulled_arm);
    }

    if (t_meta == this->meta_duration) {
      // Time to decide which algorithm has been the best, i.e. the one that has been selected
      // the most by the meta algorithm // TODO: is the selection of the best sub_alg in adapt_eve ok?
      int core_pulls = this->meta_alg->num_of_pulls[0];
      int other_pulls = this->meta_alg->num_of_pulls[1];

      if (other_pulls > core_pulls) {
        // Swap core_alg and other_alg
        MABAlgorithm *temp = this->sub_alg;
        this->sub_alg = this->other_sub_alg;
        this->other_sub_alg = temp;

        //cout << "core pulls: " << core_pulls << ", other pulls: " << other_pulls << ", keeping new alg" << endl;
      } else {
        //cout << "core pulls: " << core_pulls << ", other pulls: " << other_pulls << ", keeping old alg" << endl;
      }

      this->is_meta = false;

      // Restart CDTs
      for (int i = 0; i < this->num_of_arms; i++) {
        this->cdts[i]->reset(0);
      }
    }
  }
}
/*
GLR::GLR(string name, int num_of_arms, int M, int max_history, string sub_alg_line, boost::mt19937& rng, MAB* mab) : MetaAlgorithm(name, num_of_arms) {
  this->sub_alg = get_algorithm(sub_alg_line, &rng, mab);
  this->max_history = max_history;
  this->M = M;
  this->reset();
}

void GLR::reset(int action) {
  this->MABAlgorithm::reset(action);
  this->sub_alg->reset(action);
  if (action == -1) {
    this->rewards.clear();
    this->changepoints.clear();
    this->rewards.resize(this->num_of_arms);
    this->changepoints.assign(this->num_of_arms, 0);
    this->past_k.clear();
    this->past_k.resize(this->num_of_arms);
    for (int arm = 0; arm < this->num_of_arms; arm++) {
      this->past_k[arm].assign(1, 0);
    }
  } else {
    this->rewards[action].clear();
    this->changepoints[action] = 0;
    this->past_k[action].clear();
  }
}

int GLR::choose_action() {
  return this->sub_alg->choose_action();
}

void GLR::receive_reward(double reward, int pulled_arm) {
  this->cycle_reward(reward, pulled_arm);
  this->find_changepoints(pulled_arm);
  this->update_sub_alg(pulled_arm);
}

void GLR::cycle_reward(double new_reward, int arm_pulled) {
  this->rewards[arm_pulled].push_back(new_reward);
  if (this->rewards[arm_pulled].size() > this->max_history) {
    // Delete first element
    this->rewards[arm_pulled].erase(this->rewards[arm_pulled].begin(), this->rewards[arm_pulled].begin()+1);
  }
}

void GLR::find_changepoints(int arm_pulled) {
  if (this->rewards[arm_pulled].size() >= 2*M) {
      double p_0 = accumulate(this->rewards[arm_pulled].begin(), this->rewards[arm_pulled].begin() + this->M, 0.) / this->M;
      double p_1 = accumulate(this->rewards[arm_pulled].end() - this->M, this->rewards[arm_pulled].end(), 0.) / this->M;

      int num_of_pulls = this->rewards[arm_pulled].size();
      double sum_x_i = accumulate(this->rewards[arm_pulled].begin() + this->M, this->rewards[arm_pulled].end(), 0.);
      double l_m = sum_x_i * log(p_1/p_0) + (num_of_pulls - this->M - sum_x_i) * log((1 - p_1)/(1 - p_0));

      double prev_l = l_m;
      double highest_l = l_m;
      int worst_k = this->M;
      int best_k = this->M;
      for (int k = this->M+1; k < num_of_pulls - this->M; k++) {
        double x_k = this->rewards[arm_pulled][k-1];
        double l_k = prev_l;
        if (x_k > 0.5) { // i.e. = 1
          l_k -= log(p_1 / p_0);
        } else {
          l_k -= log((1 - p_1)/(1 - p_0));
        }
        if (l_k > highest_l) {
          highest_l = l_k;
          best_k = k;
        }
        prev_l = l_k;
      }
      this->changepoints[arm_pulled] = best_k;
  }
}

void GLR::update_sub_alg(int arm_pulled) {
  this->sub_alg->reset(arm_pulled);
  for (int i = this->changepoints[arm_pulled] + 1; i < this->rewards[arm_pulled].size(); i++) {
    this->sub_alg->receive_reward(this->rewards[arm_pulled][i], arm_pulled);
  }
}
*/

//--

Recurrent_CD_Algorithm::Recurrent_CD_Algorithm(string name, int num_of_arms, int M, double alpha, double d, string cdt_line, string sub_alg_line, bool use_history, int max_history, bool smart_resets, boost::mt19937& rng, MAB* mab) : MetaAlgorithm(name, num_of_arms) {
  this->M = M;
  this->d = d;
  this->alpha_level = alpha;
  this->max_history = max_history;
  this->global = false; 
  this->sub_alg = get_algorithm(sub_alg_line, &rng, mab);
  for (int i = 0; i < num_of_arms; i++){
    this->cdts_up.push_back(get_cdt(cdt_line + " 1"));
    this->cdts_down.push_back(get_cdt(cdt_line + " 0"));
    if (smart_resets) {
      this->sr_cdts_up.push_back(get_cdt(cdt_line + " 1"));
      this->sr_cdts_down.push_back(get_cdt(cdt_line + " 0"));
    }
  }
  this->use_history = use_history;
  this->smart_resets = smart_resets;
  this->ucb_run = false;
  if (smart_resets) {
    this->sr_change_estimates.assign(num_of_arms, 0);
  }
  this->last_resets_up.assign(num_of_arms, false);
  this->equal_concept_found.assign(num_of_arms, false);
  this->concepts.assign(num_of_arms,{}); //Assign as much concepts' containers as the arms are
  this->last_concept.assign(num_of_arms, -1); //Keep track of the similar concepts for each arm(both for single or multi)
  this->reset(-1);
}

void Recurrent_CD_Algorithm::reset(int action) {
  this->MABAlgorithm::reset(action);
  this->sub_alg->reset(action);
  this->ucb_run = false; //Set the UCB run to false
  if (action == -1) {
    for (int i = 0; i < this->num_of_arms; i++) {
      this->cdts_up[i]->reset(0);
      this->cdts_down[i]->reset(0);
      this->concepts[i].clear(); //CLEAR ALL THE CONCEPTS
      if (this->smart_resets) {
        this->sr_cdts_up[i]->reset(0);
        this->sr_cdts_down[i]->reset(0);
      }
    }
    this->till_M.clear();
    this->till_M.assign(this->num_of_arms, this->M);
    //Clear the last concept tracking
    this->last_concept.clear();
    this->last_concept.assign(this->num_of_arms, -1);
    this->collected_rewards.clear();
    this->collected_rewards.resize(this->num_of_arms);
    if (this->smart_resets) {
      this->sr_change_estimates.assign(num_of_arms, 0);
      this->sr_collected_rewards.clear();
      this->sr_collected_rewards.resize(this->num_of_arms);
    }
    this->last_resets_up.assign(num_of_arms, false);
    //Set to false all the flags for the concepts
    this->equal_concept_found.assign(num_of_arms, false);
    this->timestep = 0;
    this->to_be_reset.clear();
  } else {
    this->cdts_up[action]->reset(0);
    this->cdts_down[action]->reset(0);
    this->till_M[action] = this->M;
    this->equal_concept_found[action] = false;
    this->collected_rewards[action].clear();
    if (this->smart_resets) {
      this->sr_cdts_up[action]->reset(0);
      this->sr_cdts_down[action]->reset(0);
      this->sr_change_estimates[action] = 0;
      this->sr_collected_rewards[action].clear();
    }
    this->last_resets_up[action] = false;
    this->to_be_reset.erase(action);
  }
}

int Recurrent_CD_Algorithm::choose_action() {
  int action = -1;
  for (int i = 0; i < this->num_of_arms; i++) {
    if (this->till_M[i] > 0) {
      action = i;
      this->till_M[i]--;
      //Set the ucb_run to false if still running the sequential selection phase...
      this->ucb_run = false;
      break;
    }
  }
  if (action == -1) {
    // Then set it to true when choosing the action with any other algorithm (random or ucb1)
    this->ucb_run = true;
    action = this->sub_alg->choose_action();
  }
  return action;
}

void Recurrent_CD_Algorithm::reset_arm(int pulled_arm) {
  if(this->global){
    this->sub_alg->reset(-1);
    this->MABAlgorithm::reset(-1);
  }else{
    this->sub_alg->reset(pulled_arm);
    this->MABAlgorithm::reset(pulled_arm);
  }

  if (this->use_history) {
    int change_estimate = 0;
    bool alarm_up = this->last_resets_up[pulled_arm];
    
    if (alarm_up) change_estimate = this->cdts_up[pulled_arm]->change_estimate;
    else change_estimate = this->cdts_down[pulled_arm]->change_estimate;
    
    int history_amount = this->collected_rewards[pulled_arm].size() - change_estimate;
    int t_start = change_estimate;
    if (history_amount > this->max_history) {
      t_start = this->collected_rewards[pulled_arm].size() - this->max_history;
    }

    //Save a concept with previous information.
    double mean = 0.0;
    for(int t = 0; t != t_start; ++t){
      mean += this->collected_rewards[pulled_arm][t];
    }
    mean /= t_start;

    //SAVING THE CONCEPT
    if(this->last_concept[pulled_arm] == -1){
      Concept new_c;
      new_c.mean = mean;
      new_c.std = sqrt(mean*(1-mean));
      new_c.n = t_start; // because the rewards we want to take are those in {0,...,t_start-1}
      this->concepts[pulled_arm].push_back(new_c);
    }else{
      int last_c = this->last_concept[pulled_arm];
      //In this case we simply need to update the concept.
      double concept_mean = this->concepts[pulled_arm][last_c].mean;
      double concept_n = this->concepts[pulled_arm][last_c].prev_n;

      double new_mean = (concept_mean*concept_n + mean*t_start)/(t_start + concept_n);
      double new_std = sqrt(new_mean*(1-new_mean));
      int new_n = t_start + concept_n;

      this->concepts[pulled_arm][last_c].mean = new_mean;
      this->concepts[pulled_arm][last_c].std = new_std;
      this->concepts[pulled_arm][last_c].n = new_n;
      this->concepts[pulled_arm][last_c].prev_n = new_n;
    }

    for (int i = t_start; i < this->collected_rewards[pulled_arm].size(); i++) {
      this->sub_alg->receive_reward(this->collected_rewards[pulled_arm][i], pulled_arm);
    }
  }

  this->ucb_run = false;

  if(this->global){
    for(int arm = 0; arm != this->num_of_arms; ++arm){
      this->last_concept[arm] = -1;
      this->equal_concept_found[arm] = false;
      this->collected_rewards[arm].clear();
      this->cdts_up[arm]->reset(0);
      this->cdts_down[arm]->reset(0);
      this->till_M[arm] = this->M;
    }
  }else{
    this->last_concept[pulled_arm] = -1;
    this->equal_concept_found[pulled_arm] = false;
    this->collected_rewards[pulled_arm].clear();
    this->cdts_up[pulled_arm]->reset(0);
    this->cdts_down[pulled_arm]->reset(0);
    this->till_M[pulled_arm] = this->M;
  }
}

void Recurrent_CD_Algorithm::receive_reward(double reward, int pulled_arm) {
  this->MABAlgorithm::receive_reward(reward, pulled_arm);
  this->sub_alg->receive_reward(reward, pulled_arm);
  this->collected_rewards[pulled_arm].push_back(reward);
  
  bool new_method = false;
  bool print_bounds = false;

  // CODE FOR WRITING DOWN THE BOUNDS FOR THE MEANS.
  if(print_bounds){
    std::fstream outB("temp/exp19Bounds_"+this->name+".txt",std::ios::app);
    int tot_pulls = accumulate(this->num_of_pulls.begin(), this->num_of_pulls.end(), 0.);
    UCB1* u = dynamic_cast<UCB1*>(this->get_sub_alg());
    for(int a = 0; a != this->num_of_arms; ++a){
      double Q = u->means[a];
      double B = sqrt((2*log(tot_pulls + 1)) / max(u->num_of_pulls[a],1));
      outB << "arm " << a << " Q " << Q << " B " << B;
      outB << endl;
    }
  }

  // SIMILARITY CHECK ---------------------------------------------------------------------
  if(!new_method){
    if(!this->concepts[pulled_arm].empty()){
      for(int i = 0; i != this->concepts[pulled_arm].size(); ++i){
        double mean = this->concepts[pulled_arm][i].mean;
        double std = this->concepts[pulled_arm][i].std;
        unsigned n = this->concepts[pulled_arm][i].n;
        UCB1* ucb = dynamic_cast<UCB1*>(this->get_sub_alg());

        if(ucb == NULL){
          cout << "ERROR, algorithm does not match!!!" << endl;
          return;
        }
        double p_val = z_tost_alternative(mean, ucb->means[pulled_arm], sqrt(ucb->means[pulled_arm]*(1 - ucb->means[pulled_arm])), std, this->num_of_pulls[pulled_arm], n, -this->d, this->d);
        double num_of_concepts = this->concepts[pulled_arm].size();
        if(p_val < (this->alpha_level) && !this->equal_concept_found[pulled_arm] && this->ucb_run){
          this->last_concept[pulled_arm] = i;
          this->equal_concept_found[pulled_arm] = true;
          int new_pulls = n + this->num_of_pulls[pulled_arm];
          double new_mean = (mean*n + ucb->means[pulled_arm]*this->num_of_pulls[pulled_arm])/(new_pulls);
          double new_std = sqrt(new_mean*(1-new_mean));
          ucb->means[pulled_arm] = new_mean; //Updated means
          ucb->num_of_pulls[pulled_arm] = new_pulls; //Updated pulls
          this->num_of_pulls[pulled_arm] = new_pulls; //Updated pulls

          //UPDATE THE CONCEPT
          this->concepts[pulled_arm][i].mean = new_mean;
          this->concepts[pulled_arm][i].std = new_std;
          this->concepts[pulled_arm][i].n = new_pulls;
          break;
        }
      }
    }
  }else{
    if(!this->concepts[pulled_arm].empty()){

      int c = -1;
      double min_pval = 1.0;
      UCB1* ucb = dynamic_cast<UCB1*>(this->get_sub_alg());

      if(ucb == NULL){
        cout << "ERROR, algorithm does not match!!!" << endl;
        return;
      }

      for(int i = 0; i != this->concepts[pulled_arm].size(); ++i){
        double mean = this->concepts[pulled_arm][i].mean;
        double std = this->concepts[pulled_arm][i].std;
        unsigned n = this->concepts[pulled_arm][i].n;
        double p_val = z_tost(ucb->means[pulled_arm], mean, std, this->num_of_pulls[pulled_arm],-this->d,this->d);
        double num_of_concepts = this->concepts[pulled_arm].size();

        if(p_val < min_pval){
          min_pval = p_val;
          c = i;
        }
      }

      if(min_pval < (this->alpha_level) && !this->equal_concept_found[pulled_arm] && this->ucb_run){
        double mean = this->concepts[pulled_arm][c].mean;
        double std = this->concepts[pulled_arm][c].std;
        unsigned n = this->concepts[pulled_arm][c].n;
        //RESTORE THE CONCEPT!
        // std::cout << "detected! " << this->timestep << " for arm " << pulled_arm <<std::endl;
        // std::cout << "mean: " << mean << std::endl;
        // std::cout << "std: " << std << std::endl;
        // std::cout << "n: " << n << std::endl;
        // std::cout << "mean(" << pulled_arm << "): "<<  ucb->means[pulled_arm] << " " << this->num_of_pulls[pulled_arm] <<  std::endl;
        this->last_concept[pulled_arm] = c;
        this->equal_concept_found[pulled_arm] = true;
        int new_pulls = n + this->num_of_pulls[pulled_arm];
        double new_mean = (mean*n + ucb->means[pulled_arm]*this->num_of_pulls[pulled_arm])/(new_pulls);
        double new_std = sqrt(new_mean*(1-new_mean));
        
        // ---- insert here the restoration part ----
        // Solutions: 
        // 1) through the collected reward i insert again all the reward acquired in the previous run
        // 2) I make something to directly modify the means in the UCB1 algorithm.  
        // 3) Directly take the UCB1 sub-algorithm with a dynamic cast
        //UCB1* ucb = dynamic_cast<UCB1*>(this->get_sub_alg());
        ucb->means[pulled_arm] = new_mean;
        ucb->num_of_pulls[pulled_arm] = new_pulls;
        this->num_of_pulls[pulled_arm] = new_pulls;

        //UPDATE THE CONCEPT
        this->concepts[pulled_arm][c].mean = new_mean;
        this->concepts[pulled_arm][c].std = new_std;
        this->concepts[pulled_arm][c].n = new_pulls;
      } 
    }
  }

  // ---------------------------------------------------------------------------

  bool alarm_up = this->cdts_up[pulled_arm]->run(reward);
  bool alarm_down = this->cdts_down[pulled_arm]->run(reward);
  bool alarm = alarm_up || alarm_down;
  if (alarm_up) this->last_resets_up[pulled_arm] = true;
  else if (alarm_down) this->last_resets_up[pulled_arm] = false;

  // CHANGE DETECTION PART
  if (alarm) {
    if (!this->smart_resets) {
      this->reset_arm(pulled_arm);
    } else {
      // Directly check UCB1...
      UCB1* u = dynamic_cast<UCB1*>(this->get_sub_alg());
      int best_action = -1;
      if (u != NULL) {
        best_action = distance(u->means.begin(), max_element(u->means.begin(), u->means.end()));
      }

      this->to_be_reset.insert(pulled_arm);

      if (pulled_arm == best_action && alarm_down) {
        for (auto a : this->to_be_reset) {
          this->reset_arm(a);
        }
        this->to_be_reset.clear();
      } else if (pulled_arm != best_action && alarm_up) {
        vector<int> arms_to_reset;
        arms_to_reset.push_back(pulled_arm);
        if (this->to_be_reset.find(best_action) != this->to_be_reset.end()) {
          arms_to_reset.push_back(best_action);
        }
        for (auto a : arms_to_reset) {
          this->reset_arm(a);
          this->to_be_reset.erase(a);
        }
      } else {
        this->cdts_up[pulled_arm]->reset(1);
        this->cdts_down[pulled_arm]->reset(1);
      }
    }
  }
  this->timestep++;
}


//// ________________________ GLOBAL ALGORITHM ________________________________


Recurrent_Global_CD_Algorithm::Recurrent_Global_CD_Algorithm(string name, int num_of_arms, int M, double alpha, double d, string cdt_line, string sub_alg_line, bool use_history, int max_history, bool smart_resets, boost::mt19937& rng, MAB* mab) : MetaAlgorithm(name, num_of_arms) {
  this->M = M;
  this->d = d;
  this->alpha_level = alpha;
  this->max_history = max_history;
  this->sub_alg = get_algorithm(sub_alg_line, &rng, mab);
  for (int i = 0; i < num_of_arms; i++){
    this->cdts_up.push_back(get_cdt(cdt_line + " 1"));
    this->cdts_down.push_back(get_cdt(cdt_line + " 0"));
    if (smart_resets) {
      this->sr_cdts_up.push_back(get_cdt(cdt_line + " 1"));
      this->sr_cdts_down.push_back(get_cdt(cdt_line + " 0"));
    }
  }
  this->use_history = use_history;
  this->smart_resets = smart_resets;
  this->ucb_run = false;
  if (smart_resets) {
    this->sr_change_estimates.assign(num_of_arms, 0);
  }
  this->last_resets_up.assign(num_of_arms, false);
  this->equal_concept_found = false;
  this->last_concept = -1;
  this->reset(-1);
}


void Recurrent_Global_CD_Algorithm::reset(int action) {
  this->MABAlgorithm::reset(action);
  this->sub_alg->reset(action);
  this->ucb_run = false; //Set the UCBS run to false
  if (action == -1) {
    for (int i = 0; i < this->num_of_arms; i++) {
      this->cdts_up[i]->reset(0);
      this->cdts_down[i]->reset(0);
      this->concepts.clear(); //CLEAR ALL THE CONCEPTS
      if (this->smart_resets) {
        this->sr_cdts_up[i]->reset(0);
        this->sr_cdts_down[i]->reset(0);
      }
    }
    this->till_M.clear();
    this->till_M.assign(this->num_of_arms, this->M);
    //Clear the last concept tracking
    this->last_concept = -1;
    //Clear collected rewards and history for arms
    this->arm_pulls_history.clear();
    this->collected_rewards.clear();
    this->collected_rewards.resize(this->num_of_arms);
    if (this->smart_resets) {
      this->sr_change_estimates.assign(num_of_arms, 0);
      this->sr_collected_rewards.clear();
      this->sr_collected_rewards.resize(this->num_of_arms);
    }
    this->last_resets_up.assign(num_of_arms, false);
    //Set to false all the flags for the concepts
    this->equal_concept_found = false;
    this->timestep = 0;
    this->to_be_reset.clear();
  } else {
    this->cdts_up[action]->reset(0);
    this->cdts_down[action]->reset(0);
    this->till_M[action] = this->M;
    this->equal_concept_found = false;
    this->collected_rewards[action].clear();
    if (this->smart_resets) {
      this->sr_cdts_up[action]->reset(0);
      this->sr_cdts_down[action]->reset(0);
      this->sr_change_estimates[action] = 0;
      this->sr_collected_rewards[action].clear();
    }
    this->last_resets_up[action] = false;
    this->to_be_reset.erase(action);
  }
}

int Recurrent_Global_CD_Algorithm::choose_action() {
  int action = -1;
  for (int i = 0; i < this->num_of_arms; i++) {
    if (this->till_M[i] > 0) {
      action = i;
      this->till_M[i]--;
      //Set the ucb_run to false if still running the round phase...
      this->ucb_run = false;
      break;
    }
  }
  if (action == -1) {
    // Then set it to true when choosing the action with any other algorithm (random or ucb1)
    this->ucb_run = true;
    action = this->sub_alg->choose_action();
  }
  //Insert the chosen action in the 
  this->arm_pulls_history.push_back(action);
  return action;
}

vector<unsigned> Recurrent_Global_CD_Algorithm::num_pulls_before(int t_start, int action, int& pos){
  vector<unsigned> pulls(this->num_of_arms, 0);
  pos = 0;
  for(auto& el : this->arm_pulls_history){
    pulls[el] +=1;
    ++pos;
    if(el==action && pulls[el]==t_start) break;
  }
  return pulls;
}

vector<unsigned> Recurrent_Global_CD_Algorithm::num_pulls_after(int pos){
  vector<unsigned> pulls(this->num_of_arms, 0);
  for(int t = pos; t != this->arm_pulls_history.size(); ++t){
    pulls[this->arm_pulls_history[t]]+=1;
  }
  return pulls;
}

void Recurrent_Global_CD_Algorithm::reset_arm(int pulled_arm) {
  //Reset everything
  this->sub_alg->reset(-1); //Resets the values for the possible subalgs.
  this->MABAlgorithm::reset(-1); //Resets the number of pulls for the arms.
  
  if (this->use_history) {
    int change_estimate = 0;
    bool alarm_up = this->last_resets_up[pulled_arm];
    
    if (alarm_up) change_estimate = this->cdts_up[pulled_arm]->change_estimate;
    else change_estimate = this->cdts_down[pulled_arm]->change_estimate;
    
    int history_amount = this->collected_rewards[pulled_arm].size() - change_estimate;
    int t_start = change_estimate;
    if (history_amount > this->max_history) {
      t_start = this->collected_rewards[pulled_arm].size() - this->max_history;
    }

    //Save a concept with previous information.
    vector<double> ms(this->num_of_arms, 0.0);
    int pos = 0;
    vector<unsigned> ps = this->num_pulls_before(t_start, pulled_arm, pos);
    vector<unsigned> ps_after = this->num_pulls_after(pos);
    for(int a = 0; a != this->num_of_arms; ++a){
      double mean = 0.0;
      for(int t = 0; t < ps[a]; ++t){
        mean += this->collected_rewards[a][t];
      }
      if(ps[a]>0){
        mean /= ps[a];
      }
      ms[a] = mean;
    }

    //SAVING THE CONCEPT
    if(this->last_concept == -1){
      multiConcept new_c(this->num_of_arms);
      new_c.mean = ms;
      new_c.std.assign(this->num_of_arms,0.0);
      new_c.n = ps;
      new_c.prev_n = ps;
      for(int i = 0; i != this->num_of_arms; ++i){
        new_c.std[i] = sqrt(new_c.mean[i]*(1 - new_c.mean[i]));
      }
      // Insert new concept into the concepts of the pulled arm
      // std::cout << "Concept saved" << std::endl;
      // new_c.print();
      this->concepts.push_back(new_c);
    }else{
      int last_c = this->last_concept;
      //In this case we simply need to update the concept.
      
      for(int a = 0; a != this->num_of_arms; ++a){
        this->concepts[last_c].mean[a] = (this->concepts[last_c].mean[a]*this->concepts[last_c].n[a] + ms[a]*ps[a])/(this->concepts[last_c].n[a] + ps[a]);
        this->concepts[last_c].std[a] = sqrt(this->concepts[last_c].mean[a]*(1 - this->concepts[last_c].mean[a]));
        this->concepts[last_c].n[a] = this->concepts[last_c].prev_n[a] + ps[a];
        this->concepts[last_c].prev_n[a] = this->concepts[last_c].n[a];
      }
      // std::cout << "Concept updated (" << last_c << ") " << std::endl;
      // this->concepts[last_c].print();
    }
    for(int a = 0; a != this->num_of_arms; ++a){
      for (int i = ps_after[a]; i < this->collected_rewards[a].size(); i++) {
        this->sub_alg->receive_reward(this->collected_rewards[a][i], a);
        //UPDATE THE MEANS WITH THE HISTORY.
      }
    }
  }

  this->ucb_run = false;
  this->equal_concept_found = false;
  this->last_concept = -1;
  this->arm_pulls_history.clear();
 
  for(int arm = 0; arm != this->num_of_arms; ++arm){
      this->collected_rewards[arm].clear();
      this->cdts_up[arm]->reset(0);
      this->cdts_down[arm]->reset(0);
      this->till_M[arm] = this->M;
  }
}


void Recurrent_Global_CD_Algorithm::receive_reward(double reward, int pulled_arm){
  this->MABAlgorithm::receive_reward(reward, pulled_arm);
  this->sub_alg->receive_reward(reward, pulled_arm);
  this->collected_rewards[pulled_arm].push_back(reward);

  // SIMILARITY CHECK ---------------------------------------------------------------------

  if(!this->concepts.empty()){
    for(int i = 0; i != this->concepts.size(); ++i){
      UCB1* ucb = dynamic_cast<UCB1*>(this->get_sub_alg());

      if(ucb == NULL){
        cout << "ERROR, algorithm does not match!!!" << endl;
        return;
      }
      
      double max_p_val = 0.0;
      for(int a = 0; a != this->num_of_arms; ++a){
        //double p_val = z_tost(ucb->means[a], this->concepts[i].mean[a], this->concepts[i].std[a], this->num_of_pulls[a],-this->d ,this->d);
        double p_val = z_tost_alternative(ucb->means[a], this->concepts[i].mean[a], sqrt(ucb->means[a]*(1-ucb->means[a])),this->concepts[i].std[a],this->num_of_pulls[a],this->concepts[i].n[a],-this->d,this->d);
        max_p_val = max(max_p_val, p_val);
      }
      double num_of_concepts = this->concepts.size();
      if(max_p_val < (this->alpha_level) && !this->equal_concept_found && this->ucb_run){
        //RESTORE THE CONCEPT!
        // std::cout << "detected! " << this->timestep << " for arm " << pulled_arm <<std::endl;
        // this->concepts[i].print();
        // std::cout << "Real means: " << std::endl;
        // for(int a = 0; a != this->num_of_arms; ++a){
        //   std::cout << "arm " << a << " -> mean " << ucb->means[a] << "; n " << this->num_of_pulls[a] << std::endl;
        // }
        this->last_concept = i;
        this->equal_concept_found = true;

        //Update the concept...
        for(int a = 0; a != this->num_of_arms; ++a){
          this->concepts[i].mean[a] = (this->concepts[i].mean[a]*this->concepts[i].n[a] + ucb->means[pulled_arm]*this->num_of_pulls[pulled_arm])/(this->num_of_pulls[pulled_arm] + this->concepts[i].n[a]);
          this->concepts[i].std[a] = sqrt(this->concepts[i].mean[a]*(1 - this->concepts[i].mean[a]));
          this->concepts[i].n[a] += this->num_of_pulls[pulled_arm];
        
        
        //... and update the means of the algorithm
          ucb->means[a] = this->concepts[i].mean[a];
          ucb->num_of_pulls[a] = this->concepts[i].n[a];
          this->num_of_pulls[a] = this->concepts[i].n[a];

        }
        break;
      }
    }
  }

  // ---------------------------------------------------------------------------

  bool alarm_up = this->cdts_up[pulled_arm]->run(reward);
  bool alarm_down = this->cdts_down[pulled_arm]->run(reward);
  bool alarm = alarm_up || alarm_down;
  if (alarm_up) this->last_resets_up[pulled_arm] = true;
  else if (alarm_down) this->last_resets_up[pulled_arm] = false;

  // CHANGE DETECTION PART
  if (alarm) {
    //std::cout << "Detection at: " << this->timestep << std::endl;
    if (!this->smart_resets) {
      this->reset_arm(pulled_arm);
    } else {
      // Directly check UCB1...
      UCB1* u = dynamic_cast<UCB1*>(this->get_sub_alg());
      int best_action = -1;
      if (u != NULL) {
        best_action = distance(u->means.begin(), max_element(u->means.begin(), u->means.end()));
      }

      this->to_be_reset.insert(pulled_arm);

      if (pulled_arm == best_action && alarm_down) {
        for (auto a : this->to_be_reset) {
          this->reset_arm(a);
        }
        this->to_be_reset.clear();
      } else if (pulled_arm != best_action && alarm_up) {
        vector<int> arms_to_reset;
        arms_to_reset.push_back(pulled_arm);
        if (this->to_be_reset.find(best_action) != this->to_be_reset.end()) {
          arms_to_reset.push_back(best_action);
        }
        for (auto a : arms_to_reset) {
          this->reset_arm(a);
          this->to_be_reset.erase(a);
        }
      } else {
        this->cdts_up[pulled_arm]->reset(1);
        this->cdts_down[pulled_arm]->reset(1);
      }
    }
  }
  this->timestep++;
}



// MONITORED_UCB

Monitored_UCB::Monitored_UCB(string name, int num_of_arms, unsigned windowSize, double thresh, string sub_alg_line, bool use_history, boost::mt19937& rng, MAB* mab):MetaAlgorithm(name,num_of_arms){
  this->windowSize = windowSize;
  for(int i = 0; i != num_of_arms; i++){
    this->monitors.emplace_back(windowSize, thresh);
  }
  this->use_history = use_history;
  this->sub_alg = get_algorithm(sub_alg_line, &rng, mab);
  this->reset(-1);
}

void Monitored_UCB::reset_arm(int pulled_arm){
  this->sub_alg->reset(pulled_arm);
  this->MABAlgorithm::reset(pulled_arm);
  
  if(this->use_history){
    //Insert Code for the managing of the story
  }

  this->monitors[pulled_arm].reset(0);
}

void Monitored_UCB::reset(int action){
  this->sub_alg->reset(action);
  this->MABAlgorithm::reset(action);
  if(action == -1){
    for(int i = 0; i != this->num_of_arms; ++i){
      this->monitors[i].reset(0);
    }
  }else{
    this->monitors[action].reset(0);
  }
}

int Monitored_UCB::choose_action(){
  return this->sub_alg->choose_action();
}

void Monitored_UCB::receive_reward(double reward, int pulled_arm){
  this->MABAlgorithm::receive_reward(reward, pulled_arm);
  this->sub_alg->receive_reward(reward, pulled_arm);

  if(this->monitors[pulled_arm].run(reward)){
    this->reset(-1);
  }
}


