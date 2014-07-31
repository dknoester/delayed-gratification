/* delay.h
 *
 * This file is part of the Himalaya project.
 *
 * Copyright 2014 David B. Knoester.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _DELAY_H_
#define _DELAY_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

#include <ea/datafile.h>
#include <ea/metadata.h>


using namespace ealib;

LIBEA_MD_DECL(DELAY_GENERATIONS, "delay.generations", int);
LIBEA_MD_DECL(DELAY_W_REAL, "delay.w_real", double);
LIBEA_MD_DECL(DELAY_W_EFF, "delay.w_eff", double);
LIBEA_MD_DECL(DELAY_RANDOM_INSERT, "delay.random_insert", double);

struct delayed_priority {
    template <typename EA>
    double operator()(typename EA::individual_type& ind, EA& ea) {
        typedef typename EA::individual_ptr_type individual_ptr_type;
        
        double w = ind.priority();
        put<DELAY_W_REAL>(w, ind);
        
        if(get<IND_GENERATION>(ind) > 0) {
            individual_ptr_type p=ind.traits().lod_parent();
            for(int i=1; (i<=get<DELAY_GENERATIONS>(ea)) && (get<IND_GENERATION>(*p) >= 0); ++i) {
                w = get<DELAY_W_REAL>(*p);
                p = p->traits().lod_parent();
            }
        }
        
        put<DELAY_W_EFF>(w, ind);
        return w;
    }
};


/*! Delay the fitness of an indivdual by up to DELAY_GENERATIONS number of
 ancestors along its lineage.
 
 For example, the effective fitness w_eff of individual i is the real fitness
 w_real of its n'th ancestor.
 */
template <typename FitnessFunction>
struct generation_delay : public FitnessFunction {
    typedef FitnessFunction parent;
    
    template <typename Individual, typename EA>
    double operator()(Individual& ind, EA& ea) {
        typedef typename EA::individual_ptr_type individual_ptr_type;
        
        double w = parent::operator()(ind,ea);
        put<DELAY_W_REAL>(w, ind);
        
        if(get<IND_GENERATION>(ind) > 0) {
            individual_ptr_type p=ind.traits().lod_parent();
            for(int i=1; (i<=get<DELAY_GENERATIONS>(ea)) && (get<IND_GENERATION>(*p) >= 0); ++i) {
                w = get<DELAY_W_REAL>(*p);
                p = p->traits().lod_parent();
            }
        }
        
        put<DELAY_W_EFF>(w, ind);
        return w;
    }
};


/*! Datafile for mean generation, and mean & max fitness.
 */
template <typename EA>
struct effective_fitness : record_statistics_event<EA> {
    effective_fitness(EA& ea) : record_statistics_event<EA>(ea), _df("effective_fitness.dat") {
        _df.add_field("update")
        .add_field("mean_w_real")
        .add_field("max_w_real")
        .add_field("mean_w_eff");
    }
    
    virtual ~effective_fitness() {
    }
    
    virtual void operator()(EA& ea) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean,tag::max> > w_real;
        accumulator_set<double, stats<tag::mean> > w_eff;
        
        for(typename EA::iterator i=ea.begin(); i!=ea.end(); ++i) {
            w_real(get<DELAY_W_REAL>(*i));
            w_eff(get<DELAY_W_EFF>(*i));
        }
        
        _df.write(ea.current_update())
        .write(mean(w_real))
        .write(max(w_real))
        .write(mean(w_eff))
        .endl();
    }
    
    datafile _df;
};


#endif
