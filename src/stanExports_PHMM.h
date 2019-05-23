// Generated by rstantools.  Do not edit by hand.

/*
    Bayeshmmcts is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Bayeshmmcts is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Bayeshmmcts.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.18.1

#include <stan/model/model_header.hpp>

namespace model_PHMM_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_PHMM");
    reader.add_event(39, 37, "end", "model_PHMM");
    return reader;
}

#include <stan_meta_header.hpp>
 class model_PHMM : public prob_grad {
private:
    int N;
    vector<int> y;
    int m;
public:
    model_PHMM(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_PHMM(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;

        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_PHMM_namespace::model_PHMM";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 4;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            current_statement_begin__ = 5;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "int", context__.to_vec(N));
            validate_non_negative_index("y", "N", N);
            y = std::vector<int>(N,int(0));
            vals_i__ = context__.vals_i("y");
            pos__ = 0;
            size_t y_limit_0__ = N;
            for (size_t i_0__ = 0; i_0__ < y_limit_0__; ++i_0__) {
                y[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "m", "int", context__.to_vec());
            m = int(0);
            vals_i__ = context__.vals_i("m");
            pos__ = 0;
            m = vals_i__[pos__++];

            // validate, data variables
            current_statement_begin__ = 4;
            check_greater_or_equal(function__,"N",N,0);
            current_statement_begin__ = 5;
            for (int k0__ = 0; k0__ < N; ++k0__) {
                check_greater_or_equal(function__,"y[k0__]",y[k0__],0);
            }
            current_statement_begin__ = 6;
            check_greater_or_equal(function__,"m",m,1);
            // initialize data variables


            // validate transformed data

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 10;
            validate_non_negative_index("Gamma", "m", m);
            validate_non_negative_index("Gamma", "m", m);
            num_params_r__ += (m - 1) * m;
            current_statement_begin__ = 11;
            validate_non_negative_index("lambda", "m", m);
            num_params_r__ += m;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_PHMM() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("Gamma")))
            throw std::runtime_error("variable Gamma missing");
        vals_r__ = context__.vals_r("Gamma");
        pos__ = 0U;
        validate_non_negative_index("Gamma", "m", m);
        validate_non_negative_index("Gamma", "m", m);
        context__.validate_dims("initialization", "Gamma", "vector_d", context__.to_vec(m,m));
        std::vector<vector_d> Gamma(m,vector_d(static_cast<Eigen::VectorXd::Index>(m)));
        for (int j1__ = 0U; j1__ < m; ++j1__)
            for (int i0__ = 0U; i0__ < m; ++i0__)
                Gamma[i0__](j1__) = vals_r__[pos__++];
        for (int i0__ = 0U; i0__ < m; ++i0__)
            try {
            writer__.simplex_unconstrain(Gamma[i0__]);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable Gamma: ") + e.what());
        }

        if (!(context__.contains_r("lambda")))
            throw std::runtime_error("variable lambda missing");
        vals_r__ = context__.vals_r("lambda");
        pos__ = 0U;
        validate_non_negative_index("lambda", "m", m);
        context__.validate_dims("initialization", "lambda", "vector_d", context__.to_vec(m));
        vector_d lambda(static_cast<Eigen::VectorXd::Index>(m));
        for (int j1__ = 0U; j1__ < m; ++j1__)
            lambda(j1__) = vals_r__[pos__++];
        try {
            writer__.positive_ordered_unconstrain(lambda);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable lambda: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        typedef T__ local_scalar_t__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);

            vector<Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> > Gamma;
            size_t dim_Gamma_0__ = m;
            Gamma.reserve(dim_Gamma_0__);
            for (size_t k_0__ = 0; k_0__ < dim_Gamma_0__; ++k_0__) {
                if (jacobian__)
                    Gamma.push_back(in__.simplex_constrain(m,lp__));
                else
                    Gamma.push_back(in__.simplex_constrain(m));
            }

            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lambda;
            (void) lambda;  // dummy to suppress unused var warning
            if (jacobian__)
                lambda = in__.positive_ordered_constrain(m,lp__);
            else
                lambda = in__.positive_ordered_constrain(m);


            // transformed parameters



            // validate transformed parameters

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // model body
            {
            current_statement_begin__ = 15;
            validate_non_negative_index("log_Gamma_tr", "m", m);
            validate_non_negative_index("log_Gamma_tr", "m", m);
            vector<Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> > log_Gamma_tr(m, (Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1> (static_cast<Eigen::VectorXd::Index>(m))));
            stan::math::initialize(log_Gamma_tr, DUMMY_VAR__);
            stan::math::fill(log_Gamma_tr,DUMMY_VAR__);
            current_statement_begin__ = 16;
            validate_non_negative_index("lp", "m", m);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lp(static_cast<Eigen::VectorXd::Index>(m));
            (void) lp;  // dummy to suppress unused var warning

            stan::math::initialize(lp, DUMMY_VAR__);
            stan::math::fill(lp,DUMMY_VAR__);
            current_statement_begin__ = 17;
            validate_non_negative_index("lp_p1", "m", m);
            Eigen::Matrix<local_scalar_t__,Eigen::Dynamic,1>  lp_p1(static_cast<Eigen::VectorXd::Index>(m));
            (void) lp_p1;  // dummy to suppress unused var warning

            stan::math::initialize(lp_p1, DUMMY_VAR__);
            stan::math::fill(lp_p1,DUMMY_VAR__);


            current_statement_begin__ = 19;
            lp_accum__.add(gamma_log<propto__>(lambda, 0.10000000000000001, 0.01));
            current_statement_begin__ = 23;
            for (int i = 1; i <= m; ++i) {
                current_statement_begin__ = 24;
                for (int j = 1; j <= m; ++j) {
                    current_statement_begin__ = 25;
                    stan::model::assign(log_Gamma_tr, 
                                stan::model::cons_list(stan::model::index_uni(j), stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list())), 
                                stan::math::log(get_base1(get_base1(Gamma,i,"Gamma",1),j,"Gamma",2)), 
                                "assigning variable log_Gamma_tr");
                }
            }
            current_statement_begin__ = 27;
            stan::math::assign(lp, rep_vector(-(stan::math::log(m)),m));
            current_statement_begin__ = 29;
            for (int i = 1; i <= N; ++i) {

                current_statement_begin__ = 30;
                for (int j = 1; j <= m; ++j) {
                    current_statement_begin__ = 31;
                    stan::model::assign(lp_p1, 
                                stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list()), 
                                (log_sum_exp(add(get_base1(log_Gamma_tr,j,"log_Gamma_tr",1),lp)) + poisson_log(get_base1(y,i,"y",1),get_base1(lambda,j,"lambda",1))), 
                                "assigning variable lp_p1");
                }
                current_statement_begin__ = 33;
                stan::math::assign(lp, lp_p1);
            }
            current_statement_begin__ = 36;
            lp_accum__.add(log_sum_exp(lp));
            }

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("Gamma");
        names__.push_back("lambda");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(m);
        dims__.push_back(m);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(m);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;

        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__,params_i__);
        static const char* function__ = "model_PHMM_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector<vector_d> Gamma;
        size_t dim_Gamma_0__ = m;
        for (size_t k_0__ = 0; k_0__ < dim_Gamma_0__; ++k_0__) {
            Gamma.push_back(in__.simplex_constrain(m));
        }
        vector_d lambda = in__.positive_ordered_constrain(m);
            for (int k_1__ = 0; k_1__ < m; ++k_1__) {
                for (int k_0__ = 0; k_0__ < m; ++k_0__) {
                vars__.push_back(Gamma[k_0__][k_1__]);
                }
            }
            for (int k_0__ = 0; k_0__ < m; ++k_0__) {
            vars__.push_back(lambda[k_0__]);
            }

        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {



            // validate transformed parameters

            // write transformed parameters
            if (include_tparams__) {
            }
            if (!include_gqs__) return;
            // declare and define generated quantities



            // validate generated quantities

            // write generated quantities
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_PHMM";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_1__ = 1; k_1__ <= m; ++k_1__) {
            for (int k_0__ = 1; k_0__ <= m; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Gamma" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_0__ = 1; k_0__ <= m; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
        }


        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_1__ = 1; k_1__ <= (m - 1); ++k_1__) {
            for (int k_0__ = 1; k_0__ <= m; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "Gamma" << '.' << k_0__ << '.' << k_1__;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        for (int k_0__ = 1; k_0__ <= m; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "lambda" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (include_tparams__) {
        }


        if (!include_gqs__) return;
    }

}; // model

}

typedef model_PHMM_namespace::model_PHMM stan_model;


#endif
