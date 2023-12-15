#include "pb_class.h"
#include "GetPot"

#include <bim_distributed_vector.h>
#include <quad_operators_3d.h>
#include <mumps_class.h>
#include <lis_class.h>
// #include <lis_class_distributed.h>

#include <cmath>
#include <cstdio>
#include <fstream>

#include <p8est.h>



void
poisson_boltzmann::create_mesh ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  
  auto comp = [] (const NS::Atom &a1, const NS::Atom &a2) -> bool { return a1.radius < a2.radius; }; 
  double maxradius = std::max_element (atoms.begin (), atoms.end (), comp)->radius; 
  
  //if (mesh_shape != 2 && mesh_shape != 3)
  if (mesh_shape < 2)  
    {
      ll = atoms.begin ()->pos[0]; rr = atoms.begin ()->pos[0]; 
      l_c[0] = atoms.begin ()->pos[0]; l_c[1] = atoms.begin ()->pos[1]; l_c[2] = atoms.begin ()->pos[2]; 
      r_c[0] = atoms.begin ()->pos[0]; r_c[1] = atoms.begin ()->pos[1]; r_c[1] = atoms.begin ()->pos[2];
      auto it = [this] (const NS::Atom &a1)
        {
          for (int kk = 0; kk < 3; ++kk)
            {
              if (a1.pos[kk] > this->r_c[kk])
                this->r_c[kk] = a1.pos[kk]; 
              else if (a1.pos[kk] < this->l_c[kk])
                this->l_c[kk] = a1.pos[kk]; 
            }
        };
      std::for_each (atoms.begin (), atoms.end (), it);
    
      for (int kk = 0; kk < 3; ++kk)
        {
          if (this->rr < this->r_c[kk])
            this->rr = this->r_c[kk];
          else if (this->ll > this->l_c[kk])
            this->ll = this->l_c[kk];
        }
    
      l_c[0] -= 6*maxradius; l_c[1] -= 6*maxradius; l_c[2] -= 6*maxradius;
      r_c[0] += 6*maxradius; r_c[1] += 6*maxradius; r_c[2] += 6*maxradius;
      ll -= 6*maxradius;
      rr += 6*maxradius;
    }
  
  if (mesh_shape == 0)
    {
      if (rank == 0)
        {
          std::cout << "x: " << ll << ", " << rr << std::endl;
          std::cout << "y: " << ll << ", " << rr << std::endl;
          std::cout << "z: " << ll << ", " << rr << "\n" << std::endl;
        }
      simple_conn_num_vertices = 8;
      simple_conn_num_trees = 1;
      simple_conn_p = new double[simple_conn_num_vertices*3];
      simple_conn_t = new p4est_topidx_t[simple_conn_num_vertices];
      
      auto tmp_p = {ll, ll, ll, rr, ll, ll, ll, rr, ll, rr, rr, ll,
                   ll, ll, rr, rr, ll, rr, ll, rr, rr, rr, rr, rr};
      auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8, 1};
      
      std::copy(tmp_p.begin(), tmp_p.end(), simple_conn_p);
      std::copy(tmp_t.begin(), tmp_t.end(), simple_conn_t); 
      for (int i = 0; i<5;i++)
        bcells.push_back(std::make_pair(0, i));                  
    }
  else if (mesh_shape == 1 || mesh_shape == 2)
    {
      if (rank == 0)
        {
          std::cout << "x: " << l_c[0] << ", " << r_c[0] << std::endl;
          std::cout << "y: " << l_c[1] << ", " << r_c[1] << std::endl;
          std::cout << "z: " << l_c[2] << ", " << r_c[2] << "\n" << std::endl;
        }
      simple_conn_num_vertices = 8;
      simple_conn_num_trees = 1;
   		simple_conn_p = new double[simple_conn_num_vertices*3];
      simple_conn_t = new p4est_topidx_t[simple_conn_num_vertices];
      for (int i = 0; i<5;i++)
        bcells.push_back(std::make_pair(0, i));
      auto tmp_p = {l_c[0], l_c[1], l_c[2], 
      							r_c[0], l_c[1], l_c[2], 
      							l_c[0], r_c[1], l_c[2], 
      							r_c[0], r_c[1], l_c[2],
                    l_c[0], l_c[1], r_c[2],
                    r_c[0], l_c[1], r_c[2], 
                    l_c[0], r_c[1], r_c[2], 
                    r_c[0], r_c[1], r_c[2]};
      auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8, 1};
      
      std::copy(tmp_p.begin(), tmp_p.end(), simple_conn_p);
      std::copy(tmp_t.begin(), tmp_t.end(), simple_conn_t);
      for (int i = 0; i<5;i++)
        bcells.push_back(std::make_pair(0, i));
    }
	else if (mesh_shape == 3)
		{    
    	if (rank == 0)
        {
          std::cout << "x: " << l_c[0] << ", " << r_c[0] << std::endl;
          std::cout << "y: " << l_c[1] << ", " << r_c[1] << std::endl;
          std::cout << "z: " << l_c[2] << ", " << r_c[2] << "\n" << std::endl;
          std::cout << "Number of trees on x: " << num_trees[0] << std::endl;
          std::cout << "Number of trees on y: " << num_trees[1] << std::endl;
          std::cout << "Number of trees on z: " << num_trees[2] << "\n" << std::endl;
        }
      double bound_x = std::abs(r_c[0]-l_c[0]);
      double bound_y = std::abs(r_c[1]-l_c[1]);
      double bound_z = std::abs(r_c[2]-l_c[2]);
      double step[3] ={ bound_x/num_trees[0],
      									bound_y/num_trees[1],
      									bound_z/num_trees[2]};
      
      double bound[3] ={bound_x, bound_y, bound_z};
      										
    	make_connectivity_3d (num_trees, step, simple_conn_p, 
    												simple_conn_num_vertices, simple_conn_t,
    												simple_conn_num_trees, bcells);
    	
    	
    	for(p4est_topidx_t i =0; i < simple_conn_num_vertices; ++i)  
    		{
		  		p4est_topidx_t j = 0;
		  		simple_conn_p[3*i + j++] += l_c[0];  
		  		simple_conn_p[3*i + j++] += l_c[1];
		  		simple_conn_p[3*i + j] += l_c[2];  
		  		
    		}
    	   		 
		}
	else
		{
			std::cout << "Error during the creation of the mesh! It will be produced a cubic mesh!"<< "\n" << std::endl;
			
			if (rank == 0)
        {
					std::cout << "x: " << ll << ", " << rr << std::endl;
          std::cout << "y: " << ll << ", " << rr << std::endl;
          std::cout << "z: " << ll << ", " << rr << "\n" << std::endl;
        }
      simple_conn_num_vertices = 8;
      simple_conn_num_trees = 1;
      simple_conn_p = new double[simple_conn_num_vertices*3];
      simple_conn_t = new p4est_topidx_t[simple_conn_num_vertices];
      
      auto tmp_p = {ll, ll, ll, rr, ll, ll, ll, rr, ll, rr, rr, ll,
                   ll, ll, rr, rr, ll, rr, ll, rr, rr, rr, rr, rr};
      auto tmp_t = {1, 2, 3, 4, 5, 6, 7, 8, 1};
      
      std::copy(tmp_p.begin(), tmp_p.end(), simple_conn_p);
      std::copy(tmp_t.begin(), tmp_t.end(), simple_conn_t); 
      for (int i = 0; i<5;i++)
        bcells.push_back(std::make_pair(0, i));  
		}
  
	tmsh.read_connectivity (simple_conn_p, simple_conn_num_vertices,
                          simple_conn_t, simple_conn_num_trees);
}

double
poisson_boltzmann::levelsetfun (double x, double y, double z)
{
  double dist = 0.0;
  for (const NS::Atom& i : atoms)
    {

      dist += std::exp (surf_param * ((std::pow (x - i.pos[0], 2) +
      					                       std::pow (y - i.pos[1], 2) +
                                  	   std::pow (z - i.pos[2], 2)) /
                                      (i.radius*i.radius) - 1.0));
                                 
      if (dist > 1.5)
        break;

    }
    
  return dist;
}

double
poisson_boltzmann::is_in_ns_surf (ray_cache_t & ray_cache, double x, double y, double z)
{ 
  int rank;
  MPI_Comm_rank (mpicomm, &rank);  
  
  crossings_t & ct = ray_cache (x, z); //check if the ray is in the cache
  
  if (!ct.init && rank != 0)
    {
      std::array<double, 2> ray = {x, z};
      ray_cache.rays.erase (ray);
      return -1.;
    }
    
  int i = 0;
  if (ct.inters.size () == 0 || y < ct.inters[i])
   return 0; //if there are no inters or y_the coord is before the first intersection, the point is outside.
    
  while (i < ct.inters.size () && y > ct.inters[i]) //go on until the inters is passed
   i++;

  return (i % 2);
}

int
poisson_boltzmann::parse_options (int argc, char **argv)
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  
  GetPot g (argc, argv);
  
  if (!g.search ("--pqrfile"))
  {
     if (rank == 0)
     {
        std::cout << "Warning: No pqr file selected, using the default one." << 
        "\nTo select one use --pqrfile option followed by the desired one." << std::endl;
     }
  }

  pqrfilename = g.next ("../../data/1CCM.pqr");
  //Check that the pqr file exists
  if (rank == 0)
     std::cout << "Selected pqr file: " << pqrfilename << std::endl;
  std::ifstream pqrfile (pqrfilename);
  if (!pqrfile)
  {
     if (rank == 0)
     {
  	std::cerr << "Cannot find the pqr file" << std::endl;
  	return 1;
     }
  }
  
  if (!g.search ("--potfile"))
  {
     if (rank == 0)
     {
        std::cout << "Warning: No pot file selected, using the default one." <<
        "\nTo select one use --potfile option followed by the desired one." << std::endl;
     }
  }

  optionsfilename = g.next ("../../data/options.pot");
  //Check that the pot file exists
  if (rank == 0)
     std::cout << "Selected pot file: " << optionsfilename << std::endl;
  std::ifstream optionsfile (optionsfilename);
  if (!optionsfile)
  {
     if (rank == 0)
     {
  	std::cerr << "Cannot find the options file" << std::endl;
  	return 1;
     }
  }
  
  //Read the options from the file
  GetPot g2 (optionsfilename.c_str ());
  
  const std::string mesh_options = "mesh/";
  maxlevel = g2 ((mesh_options + "maxlevel").c_str (),  9);
  minlevel = g2 ((mesh_options + "minlevel").c_str (),  3);
  unilevel = g2 ((mesh_options + "unilevel").c_str (),  5);
  outlevel = g2 ((mesh_options + "outlevel").c_str (),  minlevel);	
  mesh_shape = g2 ((mesh_options + "mesh_shape").c_str (),  1);
  if (mesh_shape == 2)
    {
      l_c[0] = g2 ((mesh_options + "x1").c_str (),  -128.0);
      r_c[0] = g2 ((mesh_options + "x2").c_str (),  128.0);
      l_c[1] = g2 ((mesh_options + "y1").c_str (),  -128.0);
      r_c[1] = g2 ((mesh_options + "y2").c_str (),  128.0);
      l_c[2] = g2 ((mesh_options + "z1").c_str (),  -128.0);
      r_c[2] = g2 ((mesh_options + "z2").c_str (),  128.0);
      l_cr[0] = g2 ((mesh_options + "refine_x1").c_str (),  -64.0);
      r_cr[0] = g2 ((mesh_options + "refine_x2").c_str (),  64.0);
      l_cr[1] = g2 ((mesh_options + "refine_y1").c_str (),  -64.0);
      r_cr[1] = g2 ((mesh_options + "refine_y2").c_str (),  64.0);
      l_cr[2] = g2 ((mesh_options + "refine_z1").c_str (),  -64.0);
      r_cr[2] = g2 ((mesh_options + "refine_z2").c_str (),  64.0);	  
    }
	if (mesh_shape == 3)
    {
      l_c[0] = g2 ((mesh_options + "x1").c_str (),  -128.0);
      r_c[0] = g2 ((mesh_options + "x2").c_str (),  128.0);
      l_c[1] = g2 ((mesh_options + "y1").c_str (),  -128.0);
      r_c[1] = g2 ((mesh_options + "y2").c_str (),  128.0);
      l_c[2] = g2 ((mesh_options + "z1").c_str (),  -128.0);
      r_c[2] = g2 ((mesh_options + "z2").c_str (),  128.0);
      num_trees [0] = g2 ((mesh_options + "num_trees_x").c_str (),  1);
      num_trees [1] = g2 ((mesh_options + "num_trees_y").c_str (),  1);
      num_trees [2] = g2 ((mesh_options + "num_trees_z").c_str (),  1);
    }
	
  const std::string model_options = "model/";
  linearized = g2 ((model_options + "linearized").c_str (),  1);
  bc = g2 ((model_options + "bc_type").c_str (),  1);
  ionic_strength = g2 ((model_options + "ionic_strength").c_str (),  0.145);
  e_in = g2 ((model_options + "molecular_dielectric_constant").c_str (),  2.);
  e_out = g2 ((model_options + "solvent_dielectric_constant").c_str (),  78.54);
  
  const std::string surf_options = "surface/";
  int surf_type_num = g2 ((surf_options + "surface_type").c_str (),  1);
  if (surf_type_num == 1) surf_type = NS::skin;
  else if (surf_type_num == 0) surf_type = NS::ses;
  else if (surf_type_num == 2) surf_type = NS::blobby;
  else surf_type;
  surf_param = g2 ((surf_options + "surface_parameter").c_str (),  0.45);
  stern_layer = g2 ((surf_options + "stern_layer_thickness").c_str (),  2.);
  num_threads = g2 ((surf_options + "number_of_threads").c_str (),  1);
  
  const std::string alg_options = "algorithm/";
  linear_solver_name = g2 ((alg_options + "linear_solver").c_str (),  "mumps");
  linear_solver_options = g2 ((alg_options + "solver_options").c_str (),  "-i cg -p ilu [0] -tol 1.e-12");

  const std::string out_options = "output/";
  p4estfilename = g2 ((out_options + "p4estfilename").c_str (), "poisson_boltzmann_p4est");
  markerfilename = g2 ((out_options + "markerfilename").c_str (), "poisson_boltzmann_marker_0");
  surffilename = g2 ((out_options + "surffilename").c_str (), "poisson_boltzmann_surface_0");
  
  return 0;
}

void 
poisson_boltzmann::print_options ()
{
  std::cout << "\nChosen options: " << std::endl;
  std::cout << "minlevel = " << minlevel <<  "\nmaxlevel = " << maxlevel << std::endl;
  std::cout << "unilevel = " << unilevel << std::endl;
  
  std::cout << "Linearized model = " << linearized << std::endl;
  
  if (bc == 1)
    std::cout << "Dirichlet boundary conditions" << std::endl;
  else
    std::cout << "Neumann boundary conditions" << std::endl;
    
  std::cout << "e_in = " << e_in << "\ne_out = " << e_out 
            << "\nionic_strenght = " << ionic_strength << std::endl;
            
  std::cout << "Linear solver = " << linear_solver_name << std::endl;
  std::cout << "Linear solver options = " << linear_solver_options << std::endl;
  
  std::cout << "Chosen surface: " << surf_type << std::endl;
  std::cout << "Surface parameter: " << surf_param << std::endl;
  std::cout << "Stern layer thickness: " << stern_layer << std::endl;
  std::cout << "Number of threads for nanoshaper: " << num_threads << std::endl;
  if (mesh_shape == 1)
    {
      std::cout << "Mesh shape = stretched" << std::endl;
      std::cout << "with the following domain vertices: " << std::endl;
    }
  else if (mesh_shape == 0)
    {
      std::cout << "Mesh shape = cubic" << std::endl;
      std::cout << "with the following domain vertices: " << std::endl;
    }
  else 
    std::cout << "Manual setting of the mesh vertices:" << std::endl;
   
}

void
poisson_boltzmann::read_atoms_from_pqr (std::basic_istream<char> &inputfile) 
{
  static NS::Atom a;
  atoms.clear ();
  while (inputfile >> a)          
    atoms.push_back (a); 
}

void
poisson_boltzmann::write_atoms_to_pqr (std::basic_ostream<char> &outputfile) 
{
  int Atom_number = 1;

  outputfile << std::setw(10) << std::left << "fieldname" << std::setw(12) 
             << std::left <<"Atom_number" << std::setw(12) << std::left << "Atom_name" << std::setw(16) << std::left
             << "Residue_name" << std::setw(16) << std::left << "Residue_number" << std::setw(10) << std::left << "X" 
             << std::setw(10) << std::left << "Y" << std::setw(10) << std::left
             << "Z" << std::setw(10) << std::left << "Charge" << std::setw(10) << std::left << "Radius" << std::endl;

  for (auto & ii : atoms) 
    outputfile << std::setw(10) << std::left << "ATOM" << std::setw(12) 
               << std::left << Atom_number++ << std::setw(12) << std::left << ii.ai.name << std::setw(16) << std::left
               << ii.ai.resName << std::setw(16) << std::left << ii.ai.resNum << std::setw(10) << std::left << ii.pos[0] 
               << std::setw(10) << std::left << ii.pos[1] << std::setw(10) << std::left
               << ii.pos[2] << std::setw(10) << std::left << ii.charge << std::setw(10) << std::left << ii.radius << std::endl;

}

std::basic_istream<char>& 
operator>>(std::basic_istream<char>& inputfile, NS::Atom &a)
{

  int Atom_number;
  std::string Field_name;

  inputfile >> Field_name
            >> Atom_number
            >> a.ai.name
            >> a.ai.resName
            >> a.ai.resNum
            >> a.pos[0]
            >> a.pos[1]
            >> a.pos[2]
            >> a.charge
            >> a.radius;

  a.radius2 = a.radius*a.radius;
  return inputfile;
}

void
poisson_boltzmann::init_tmesh ()
{
  for (auto i = 0; i < unilevel; ++i)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (0, 1);
    }
}

void
poisson_boltzmann::init_tmesh_with_refine_box ()
{
  for (auto i = 0; i < outlevel; ++i)
    {
      tmsh.set_refine_marker (uniform_refinement);
      tmsh.refine (0, 1);
    }
	
   auto refinement = [this]
        (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 0;
            
            if (currentlevel >= this->unilevel)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {
                    if (! q->is_hanging (ii))
                      {
                        if ((q -> p(0, ii) > this->l_cr[0]) && (q -> p(0, ii) < this->r_cr[0])
			                       && (q -> p(1, ii) > this->l_cr[1]) && (q -> p(1, ii) < this->r_cr[1])
			                       && (q -> p(2, ii) > this->l_cr[2]) && (q -> p(2, ii) < this->r_cr[2]))
			                         {
			                            retval = 1;
		                              break;
			                         }
                      }
                  }
                }
                return (retval);
              };	
	
    for (auto i = 0; i < unilevel; ++i)
      {
        tmsh.set_refine_marker (refinement);
        tmsh.refine (0, 1);
      }
      
}

bool
poisson_boltzmann::is_in (const NS::Atom& i,
                          tmesh_3d::quadrant_iterator q) 
{ 
  double tol =  p4esttol * (rr-ll);
  if (mesh_shape == 2)
    tol = p4esttol * (r_c[0]-l_c[0]);
  bool retval = false;
  double l, r, t, b, f, bk;

  l = q->p (0, 0);
  r = q->p (0, 0);

  f  = q->p (1, 0);
  bk = q->p (1, 0);

  b = q->p (2, 0);
  t = q->p (2, 0);


  for (int ii = 1; ii < 8; ++ii)
    {
      l  = q->p (0, ii)  < l  ? q->p (0, ii) : l;
      r  = q->p (0, ii)  > r  ? q->p (0, ii) : r;
      f  = q->p (1, ii)  < f  ? q->p (1, ii) : f;
      bk = q->p (1, ii)  > bk ? q->p (1, ii) : bk;
      b  = q->p (2, ii)  < b  ? q->p (2, ii) : b;
      t  = q->p (2, ii)  > t  ? q->p (2, ii) : t;

    }

  retval =           (i.pos[0] > l - tol) && (i.pos[0] <= r  - tol); //make sure that the charge is assigned only once
  retval = retval && (i.pos[1] > f - tol) && (i.pos[1] <= bk - tol);
  retval = retval && (i.pos[2] > b - tol) && (i.pos[2] <= t  - tol);

  return retval;
}

void
poisson_boltzmann::refine_surface (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);
  
  int num_cycles = 2;
  if (size == 0 || surf_type == 2)
    num_cycles = 1;
  
  int coars_ref_cycles = (maxlevel - unilevel) > (unilevel - minlevel) ? (maxlevel - unilevel) : (unilevel - minlevel);
  for (int kk = 0; kk < coars_ref_cycles; ++kk)
    {
      // REFINEMENT
      {
	distributed_vector rcoeff (tmsh.num_owned_nodes ());	
		
	for (int jj = 0; jj < num_cycles; jj++)
	  {        
            ray_cache.num_req_rays = 0; //zero at each ref/coars cycle
            ray_cache.rays_list.clear (); 

            if (rank == 0 && jj == 0)
              std::cout << "Refinement: " << kk << std::endl;

            for (auto quadrant = tmsh.begin_quadrant_sweep ();
                 quadrant != tmsh.end_quadrant_sweep ();
                 ++quadrant)
              {

                for (int ii = 0; ii < 8; ++ii)
                  {
                            
                    if (! quadrant->is_hanging (ii))
                      {
                        if (surf_type == 2)
                          rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                     				                       quadrant->p (1, ii),
                                   				                         quadrant->p (2, ii));
                        else
                          {
                            rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                  					                                  quadrant->p (0, ii),
                  		            		                        quadrant->p (1, ii),
                  	         			                            quadrant->p (2, ii));
                    
                           if (rcoeff[quadrant->gt (ii)] < -0.5)
                             {
                               ray_cache.num_req_rays++;
      			        
      			                   std::array<double, 2> ray;
                               ray = {quadrant->p (0, ii), quadrant->p (2, ii)};
                               ray_cache.rays_list.insert(ray);
      			                 }
                         }
                       }
                  else
                    for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                      rcoeff[quadrant->gparent (jj, ii)] += 0.;
                  }
              }
            MPI_Barrier(mpicomm);  
            ray_cache.fill_cache();
          }
            
        auto refinement = [&rcoeff,this]
          (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 0;
            double min = 2.0; 
            double max = 0.0;
            double tmp = 0.0;

            if (currentlevel >= this->maxlevel)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {

                    if (! q->is_hanging (ii))
                      {
                        tmp = rcoeff[q->gt (ii)];

                        if (tmp > max) max = tmp;
                        if (tmp < min) min = tmp;
                      }

                  }
                if (this->surf_type == 2)
                  {
                    if (max > 1 && min < 1)
                      retval = this->maxlevel - currentlevel;
                    else
                      for (const NS::Atom& i : atoms) 
                        if (is_in (i, q))
                          {
                            retval = this->maxlevel - currentlevel;
                            break;
                          }
                   }
                else 
                  {
                    if (max > 0.5 && min < 0.5)
                      retval = this->maxlevel - currentlevel;
                    else
                      for (const NS::Atom& i : atoms) 
                        if (is_in (i, q))
                          {
                            retval = this->maxlevel - currentlevel;
                            break;
                          }
                  }
              }
	    
            return (retval);
          };

        tmsh.set_refine_marker (refinement);
        tmsh.refine (0, 1);
      }
      
      // COARSENING
      {
        distributed_vector rcoeff (tmsh.num_owned_nodes ());
        
        for (int jj = 0; jj < num_cycles; jj++)
	  {
            ray_cache.num_req_rays = 0; //zero at each ref/coars cycle
            ray_cache.rays_list.clear (); 

            if (rank == 0 && jj == 0)
              std::cout << "Coarsening: " << kk << std::endl;

            for (auto quadrant = tmsh.begin_quadrant_sweep ();
                 quadrant != tmsh.end_quadrant_sweep ();
                 ++quadrant)
              {

                for (int ii = 0; ii < 8; ++ii)
                  {
                            
                    if (! quadrant->is_hanging (ii))
                      {
                        if (surf_type == 2)
                          rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                     				 quadrant->p (1, ii),
                                   				 quadrant->p (2, ii));
                       else
                         {
                           rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                  					    quadrant->p (0, ii),
                  		            		    quadrant->p (1, ii),
                  	         			    quadrant->p (2, ii));
                    
                           if (rcoeff[quadrant->gt (ii)] < -0.5)
                             {
                               ray_cache.num_req_rays++;
      			        
      			        std::array<double, 2> ray;
                               ray = {quadrant->p (0, ii), quadrant->p (2, ii)};
      			        ray_cache.rays_list.insert(ray);
      			      }
                         }
                       }
                  else
                    for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                      rcoeff[quadrant->gparent (jj, ii)] += 0.;
                  }
              }
            MPI_Barrier(mpicomm);   
            ray_cache.fill_cache();
          }
          
        auto coarsening = [&rcoeff,this]
          (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 0;
            double min = 2.0;
            double max = 0.0;
            double tmp = 0.0;

            if (currentlevel <= this->minlevel)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {

                    if (! q->is_hanging (ii))
                      {
                        tmp = rcoeff[q->gt (ii)];

                        if (tmp > max) max = tmp;
                        if (tmp < min) min = tmp;
                      }

                   }
                 if (this->surf_type == 2)
                   {
                     if (min > 1 || max < 1)
                       retval = currentlevel - this->minlevel;  
                   }
                 else
                   {
                     if (min > 0.5 || max < 0.5)
                       retval = currentlevel - this->minlevel;
                   }

                 for (const NS::Atom& i : atoms) 
                   if (is_in (i, q))
                     {
                       retval = 0;
                       break;
                     }
               }
               
             return (retval);
          };

        tmsh.set_coarsen_marker (coarsening);
        tmsh.coarsen (0, 1);
      }
    }
}


void
poisson_boltzmann::refine_only_surface (ray_cache_t & ray_cache)
{
  int rank, size;
  MPI_Comm_size (mpicomm, &size);
  MPI_Comm_rank (mpicomm, &rank);
  
  
  int num_cycles = 2;
  if (size == 0 || surf_type == 2)
    num_cycles = 1;
  maxlevel1 = maxlevel + 4;
  int coars_ref_cycles = (maxlevel1 - unilevel) > (unilevel - minlevel1) ? (maxlevel1 - unilevel) : (unilevel - minlevel1);
  for (int kk = 0; kk < coars_ref_cycles; ++kk)
    {
      
      // REFINEMENT
      {
        distributed_vector rcoeff (tmsh.num_owned_nodes ());  
          
        for (int jj = 0; jj < num_cycles; jj++)
          {        
            ray_cache.num_req_rays = 0; //zero at each ref/coars cycle
            ray_cache.rays_list.clear (); 

            if (rank == 0 && jj == 0)
              std::cout << "Refinement: " << kk << std::endl;

            for (auto quadrant = tmsh.begin_quadrant_sweep ();
                 quadrant != tmsh.end_quadrant_sweep ();
                 ++quadrant)
              {

                for (int ii = 0; ii < 8; ++ii)
                  {
                            
                    if (! quadrant->is_hanging (ii))
                      {
                        if (surf_type == 2)
                          rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                                                   quadrant->p (1, ii),
                                                                   quadrant->p (2, ii));
                        else
                          {
                            rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                                                              quadrant->p (0, ii),
                                                              quadrant->p (1, ii),
                                                              quadrant->p (2, ii));
                    
                           if (rcoeff[quadrant->gt (ii)] < -0.5)
                             {
                               ray_cache.num_req_rays++;
                    
                               std::array<double, 2> ray;
                               ray = {quadrant->p (0, ii), quadrant->p (2, ii)};
                               ray_cache.rays_list.insert(ray);
                             }
                         }
                       }
                  else
                    for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                      rcoeff[quadrant->gparent (jj, ii)] += 0.;
                  }
              }
            MPI_Barrier(mpicomm);  
            ray_cache.fill_cache();
          }
            
        auto refinement = [&rcoeff,this]
          (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 0;
            double min = 2.0; 
            double max = 0.0;
            double tmp = 0.0;

            if (currentlevel >= this->maxlevel1)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {

                    if (! q->is_hanging (ii))
                      {
                        tmp = rcoeff[q->gt (ii)];

                        if (tmp > max) max = tmp;
                        if (tmp < min) min = tmp;
                      }

                  }
                if (this->surf_type == 2)
                  {
                    if (max > 1 && min < 1)
                      retval = this->maxlevel1 - currentlevel;
                    else
                      for (const NS::Atom& i : atoms) 
                        if (is_in (i, q))
                          {
                            retval = this->maxlevel1 - currentlevel;
                            break;
                          }
                   }
                else 
                  {
                    if (max > 0.5 && min < 0.5)
                      retval = this->maxlevel1 - currentlevel;
                    else
                      for (const NS::Atom& i : atoms) 
                        if (is_in (i, q))
                          {
                            retval = this->maxlevel1 - currentlevel;
                            break;
                          }
                  }
              }
      
            return (retval);
          };
        tmsh.set_refine_marker (refinement);
        tmsh.refine (0, 1);
      }
      

      // COARSENING
      {
        distributed_vector rcoeff (tmsh.num_owned_nodes ());

        for (int jj = 0; jj < num_cycles; jj++)
          {
            ray_cache.num_req_rays = 0; //zero at each ref/coars cycle
            ray_cache.rays_list.clear (); 

            if (rank == 0 && jj == 0)
              std::cout << "Coarsening: " << kk << std::endl;

            for (auto quadrant = tmsh.begin_quadrant_sweep ();
                 quadrant != tmsh.end_quadrant_sweep ();
                 ++quadrant)
              {

                for (int ii = 0; ii < 8; ++ii)
                  {
                            
                    if (! quadrant->is_hanging (ii))
                      {
                        if (surf_type == 2)
                          rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                                                    quadrant->p (1, ii),
                                                                   quadrant->p (2, ii));
                       else
                         {
                           rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                                                                      quadrant->p (0, ii),
                                                                      quadrant->p (1, ii),
                                                                      quadrant->p (2, ii));
                    
                           if (rcoeff[quadrant->gt (ii)] < -0.5)
                             {
                                ray_cache.num_req_rays++;
                    
                                std::array<double, 2> ray;
                                ray = {quadrant->p (0, ii), quadrant->p (2, ii)};
                                ray_cache.rays_list.insert(ray);
                             }
                         }
                       }
                  else
                    for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                      rcoeff[quadrant->gparent (jj, ii)] += 0.;
                  }
              }
            MPI_Barrier(mpicomm);   
            ray_cache.fill_cache();
          }
          
        auto coarsening = [&rcoeff,this]
          (tmesh_3d::quadrant_iterator q) -> int
          {
            int currentlevel = static_cast<int> (q->the_quadrant->level);
            int retval = 0;
            double min = 2.0;
            double max = 0.0;
            double tmp = 0.0;

            if (currentlevel <= this->minlevel1)
              retval = 0;
            else
              {
                for (int ii = 0; ii < 8; ++ii)
                  {

                    if (! q->is_hanging (ii))
                      {
                        tmp = rcoeff[q->gt (ii)];

                        if (tmp > max) max = tmp;
                        if (tmp < min) min = tmp;
                      }

                   }
                 if (this->surf_type == 2)
                   {
                     if (min > 1 || max < 1)
                       retval = currentlevel - this->minlevel1;  
                   }
                 else
                   {
                     if (min > 0.5 || max < 0.5)
                       retval = currentlevel - this->minlevel1;
                   }

                 for (const NS::Atom& i : atoms) 
                   if (is_in (i, q))
                     {
                       retval = 0;
                       break;
                     }
               }
               
             return (retval);
          };

        tmsh.set_coarsen_marker (coarsening);
        tmsh.coarsen (0, 1);
      }

    }
}



void
poisson_boltzmann::create_markers (ray_cache_t & ray_cache)
{ 

  int size, rank;
  MPI_Comm_size (mpicomm, &size); 
  MPI_Comm_rank (mpicomm, &rank); 
  
  this->marker.assign (this->tmsh.num_local_quadrants (), 0.0); //marker = 0 -> in
  
  int num_cycles = 2;
  if (size == 0 || surf_type == 2)
    num_cycles = 1;
    
  for (int jj = 0; jj < num_cycles; jj++)
    {
      ray_cache.num_req_rays = 0; //zero at each ref/coars cycle
      ray_cache.rays_list.clear (); 
      for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
           quadrant != this->tmsh.end_quadrant_sweep ();
           ++quadrant)
        {            
            int num_int_nodes = 0;
            int num_hanging = 0;
      
            for (int ii = 0; ii < 8; ++ii)
              {         
                if (! quadrant->is_hanging (ii)) 
                  {
                    if (surf_type == 2)
                      {
                        if (this->levelsetfun (quadrant->p (0, ii),
                                               quadrant->p (1, ii),
                                               quadrant->p (2, ii)) > 1.0)
                          ++num_int_nodes;
                      }
                    else
                      {
                        if (this->is_in_ns_surf (ray_cache,
                       	                quadrant->p (0, ii),
                                        quadrant->p (1, ii),
                                        quadrant->p (2, ii)) > 0.5) //inside the molecule
                          ++num_int_nodes;
                        else if (this->is_in_ns_surf (ray_cache,
                       	                     quadrant->p (0, ii),
                                             quadrant->p (1, ii),
                                             quadrant->p (2, ii)) < -0.5 )
                          {
                            ray_cache.num_req_rays++;
                            
      			                std::array<double, 2> ray;
      			                ray = {quadrant->p (0, ii), quadrant->p (2, ii)};
      			                ray_cache.rays_list.insert(ray);
      			              }
                      }
                    }
                  else
                    ++num_hanging; 
                }
              
              if (jj != 0 || num_cycles == 1)
                {
                  if (num_int_nodes == 0)  //if there's no node inside the molecule
                    this->marker[quadrant->get_forest_quad_idx ()] = 1.0; //quadrant is out 
                  else if (num_int_nodes < (8 - num_hanging)) //if the non hanging nodes are not all inside
                    this->marker[quadrant->get_forest_quad_idx ()] = 1.0/2.0; //"border"
                  //else: all the nodes are inside: the quadrant is inside and the marker value is 0
                }
            }
        MPI_Barrier(mpicomm);     
        ray_cache.fill_cache();
      }          
}

void
poisson_boltzmann::export_tmesh (ray_cache_t & ray_cache)
{
    int size, rank;
    MPI_Comm_size (mpicomm, &size);
    MPI_Comm_rank (mpicomm, &rank);
    
    distributed_vector rcoeff (tmsh.num_owned_nodes ());
    
    for (auto quadrant = tmsh.begin_quadrant_sweep ();
         quadrant != tmsh.end_quadrant_sweep ();
         ++quadrant)
      {

        for (int ii = 0; ii < 8; ++ii)
          {
                              
            if (! quadrant->is_hanging (ii))
              {
                if (surf_type == 2)
                  rcoeff[quadrant->gt (ii)] = levelsetfun (quadrant->p (0, ii),
                                                           quadrant->p (1, ii),
                                                           quadrant->p (2, ii));
                else
                  {
                    rcoeff[quadrant->gt (ii)] = is_in_ns_surf (ray_cache,
                                                               quadrant->p (0, ii),
                                                               quadrant->p (1, ii),
                                                               quadrant->p (2, ii));
                  }
               }
                      
             else
               for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                 rcoeff[quadrant->gparent (jj, ii)] += 0.;
           }
      }
    
    bim3a_solution_with_ghosts (tmsh, rcoeff, replace_op);
    tmsh.octbin_export (surffilename.c_str (), rcoeff);
}

void
poisson_boltzmann::export_marked_tmesh ()
{
  tmsh.octbin_export_quadrant (markerfilename.c_str (), marker);
}

void
poisson_boltzmann::export_p4est ()
{
  tmsh.save (p4estfilename.c_str ());
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

// 4-points Gauss quadature nodes and weights (in [0, 1]).
static constexpr double gn[4] =
  {6.94318442029737e-02, 3.30009478207572e-01,
   6.69990521792428e-01, 9.30568155797026e-01};
static constexpr double gw[4] =
  {1.73927422568727e-01, 3.26072577431273e-01,
   3.26072577431273e-01, 1.73927422568727e-01};

// Transform nodes from [0, 1] to [x[0], x[1]].
static inline double
xformx (const double *x, const double X)
{ return (X * (x[1] - x[0]) + x[0]); }

// Transform weights from [0, 1] to [x[0], x[1]].
static inline double
xformw (const double *x, const double w)
{ return (w * (x[1] - x[0])); }


// Approximate integral of fun on [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]].
static double
quad_integral (const double *x, const double *y, const double *z,
               std::function<double (double, double, double)> fun)
{
  int ix, jy, kz;
  double wx = 0, nx = 0, sum = 0;
  for (ix = 0; ix < 4; ++ix)
    {
      wx = xformw (x, gw[ix]);
      nx = xformx (x, gn[ix]);
      double wy = 0, ny = 0;
      for (jy = 0; jy < 4; ++jy)
        {
          wy = xformw (y, gw[jy]);
          ny = xformx (y, gn[jy]);
          for (kz = 0; kz < 4; ++kz)
            sum += fun (nx, ny, xformx(z, gn[kz])) *
              wx * wy * xformw(z, gw[kz]);
        }
    }
  return (sum);
}


// Evaluate u (using Q1 basis functions
// on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
q1 (double X, double Y, double Z, const double *x,
    const double *y, const double *z, const double *u)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);
  double Xx1 = (X - x[1]);
  double Xx0 = (X - x[0]);
  double Yy1 = (Y - y[1]);
  double Yy0 = (Y - y[0]);
  double Zz1 = (Z - z[1]);
  double Zz0 = (Z - z[0]);

  double num = u[0] * -Xx1  * Yy1 * Zz1 +
    u[1] * Xx0 * Yy1 * Zz1 +
    u[2] * Xx1 * Yy0 * Zz1 +
    u[3] * -Xx0  * Yy0 * Zz1 +
    u[4] * Xx1  * Yy1 * Zz0 +
    u[5] * -Xx0 * Yy1 * Zz0 +
    u[6] * -Xx1 * Yy0 * Zz0 +
    u[7] * Xx0  * Yy0 * Zz0;

  return (num / hxhyhz);
}


// Evaluate Nedelec x-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
dudx (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double db1 = (u[1] - u[0]);
  double dt1 = (u[3] - u[2]);

  double db2 = (u[5] - u[4]);
  double dt2 = (u[7] - u[6]);

  double result = (db1 * (y[1] - Y) * (z[1] - Z) +
                   dt1 * (Y - y[0]) * (z[1] - Z) +
                   db2 * (y[1] - Y) * (Z - z[0]) +
                   dt2 * (Y - y[0]) * (Z - z[0])) / hxhyhz;

  return  result;
}

// Evaluate x-component of dielectric displacement
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
Dx (double X, double Y, double Z, const double *x,
    const double *y, const double *z, const double *u, const double *eps)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double db1 = (u[1] - u[0])*2.0*(1.0/eps[1] + 1.0/eps[0]);
  double dt1 = (u[3] - u[2])*2.0*(1.0/eps[3] + 1.0/eps[2]);

  double db2 = (u[5] - u[4])*2.0*(1.0/eps[5] + 1.0/eps[4]);
  double dt2 = (u[7] - u[6])*2.0*(1.0/eps[7] + 1.0/eps[6]);

  double result = (db1 * (y[1] - Y) * (z[1] - Z) +
                   dt1 * (Y - y[0]) * (z[1] - Z) +
                   db2 * (y[1] - Y) * (Z - z[0]) +
                   dt2 * (Y - y[0]) * (Z - z[0])) / hxhyhz;

  return  result;
}

// Evaluate Nedelec y-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
dudy (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double dl1 = (u[2] - u[0]);
  double dr1 = (u[3] - u[1]);

  double dl2 = (u[6] - u[4]);
  double dr2 = (u[7] - u[5]);

  double result = (dl1 * (x[1] - X) * (z[1] - Z) +
                   dr1 * (X - x[0]) * (z[1] - Z) +
                   dl2 * (x[1] - X) * (Z - z[0]) +
                   dr2 * (X - x[0]) * (Z - z[0])) / hxhyhz;

  return result;
}

// Evaluate y-component of dielectric displacement
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
Dy (double X, double Y, double Z, const double *x,
    const double *y, const double *z, const double *u, const double *eps)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double dl1 = (u[2] - u[0])*2.0*(1.0/eps[2] + 1.0/eps[0]);
  double dr1 = (u[3] - u[1])*2.0*(1.0/eps[3] + 1.0/eps[1]);

  double dl2 = (u[6] - u[4])*2.0*(1.0/eps[6] + 1.0/eps[4]);
  double dr2 = (u[7] - u[5])*2.0*(1.0/eps[7] + 1.0/eps[5]);

  double result = (dl1 * (x[1] - X) * (z[1] - Z) +
                   dr1 * (X - x[0]) * (z[1] - Z) +
                   dl2 * (x[1] - X) * (Z - z[0]) +
                   dr2 * (X - x[0]) * (Z - z[0])) / hxhyhz;

  return result;
}

// Evaluate Nedelec z-gradient of u
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
dudz (double X, double Y, double Z, const double *x,
      const double *y, const double *z, const double *u)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double d11 = (u[4] - u[0]);
  double d21 = (u[5] - u[1]);

  double d12 = (u[6] - u[2]);
  double d22 = (u[7] - u[3]);

  double result = (d11 * (x[1] - X) * (y[1] - Y) +
                   d21 * (X - x[0]) * (y[1] - Y) +
                   d12 * (x[1] - X) * (Y - y[0]) +
                   d22 * (X - x[0]) * (Y - y[0])) / hxhyhz;

  return result;
}

// Evaluate z-component of dielectric displacement
// (on quadrant [x[0], x[1]] x [y[0], y[1]] x [z[0], z[1]])
// at (X, Y, Z).
static double
Dz (double X, double Y, double Z, const double *x,
    const double *y, const double *z, const double *u, const double *eps)
{
  double hxhyhz = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0]);

  double d11 = (u[4] - u[0])*2.0*(1.0/eps[4] + 1.0/eps[0]);
  double d21 = (u[5] - u[1])*2.0*(1.0/eps[5] + 1.0/eps[1]);

  double d12 = (u[6] - u[2])*2.0*(1.0/eps[6] + 1.0/eps[2]);
  double d22 = (u[7] - u[3])*2.0*(1.0/eps[7] + 1.0/eps[3]);

  double result = (d11 * (x[1] - X) * (y[1] - Y) +
                   d21 * (X - x[0]) * (y[1] - Y) +
                   d12 * (x[1] - X) * (Y - y[0]) +
                   d22 * (X - x[0]) * (Y - y[0])) / hxhyhz;

  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



/*

void
poisson_boltzmann::mumps_compute_electric_potential ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  if (rank == 0)
     std::cout << "\nStarting MUMPS solution" << std::endl;
     
  // diffusion
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/(e*e);   //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/(e*e); //adim e_out
  epsilon.assign (tmsh.num_local_quadrants (), eps_in); //e_in
  epsilon_in.assign (tmsh.num_local_quadrants (), 0.0); 
  epsilon_out = epsilon_in;
  
  distributed_vector  psi (tmsh.num_owned_nodes ());
  psi.get_owned_data ().assign (psi.get_owned_data ().size (), 0.0);
  
  std::vector<double> ones_in (tmsh.num_local_quadrants (), 0.0);
  
  for (auto epsp = epsilon.begin (), mp = marker.begin ();
       epsp != epsilon.end () || mp != marker.end ();
       ++epsp, ++mp) 
    if ((*mp) < 0.6)
      (*epsp) = eps_in; 
    else  
      (*epsp) = eps_out;

  for (auto epsp_in = epsilon_in.begin (), 
            epsp_out = epsilon_out.begin (), 
            onesp_in = ones_in.begin (), 
            mp = marker.begin ();
       epsp_in != epsilon_in.end () 
       || epsp_out != epsilon_out.end () 
       || onesp_in != ones_in.end ()
       || mp != marker.end ();
       ++epsp_in, ++epsp_out, ++onesp_in, ++mp) 
    if ((*mp) < 0.6){
      (*epsp_in) = eps_in;
      (*onesp_in) = 1.0; 
      }
    else  
      (*epsp_out) = eps_out;
  
  //reaction 
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/(e_0*e_out*kb*T);  
  reaction.assign (tmsh.num_local_quadrants (), 0.0);
  
  distributed_vector ones (tmsh.num_owned_nodes ()); 
  ones.get_owned_data ().assign (ones.get_owned_data ().size (), 1.0); 
  
  for (auto rp = reaction.begin (), mp = marker.begin ();
       rp != reaction.end () || mp != marker.end ();
       ++rp, ++mp)
    if ((*mp) > 0.6) 
      (*rp) = eps_out*k2; 

  tmsh.octbin_export_quadrant ("epsilon_0", epsilon);
  tmsh.octbin_export_quadrant ("reaction_0", reaction);
  
  distributed_vector rho_fixed (tmsh.num_owned_nodes ()); 
  rho_fixed.get_owned_data ().assign (rho_fixed.get_owned_data ().size (), 0.0); 
  std::vector<double> const_ones (tmsh.num_local_quadrants (), 1.0);
  
  distributed_vector  vol_patch (tmsh.num_owned_nodes (), mpicomm); 
  bim3a_solution_with_ghosts (tmsh, ones, replace_op);
  bim3a_rhs (tmsh, const_ones, ones, vol_patch);
  vol_patch.assemble ();
  
  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
           quadrant != this->tmsh.end_quadrant_sweep ();
           ++quadrant)
    {
      for (const NS::Atom& i : atoms) 
        if (is_in (i, quadrant)) 
          {
            //linear approx:
            double volume = (quadrant->p(0, 7) - quadrant->p(0, 0)) *
                            (quadrant->p(1, 7) - quadrant->p(1, 0)) *
                            (quadrant->p(2, 7) - quadrant->p(2, 0)); //volume
            for (int ii = 0; ii < 8; ++ii)
              {
                double weigth = std::abs ((i.pos[0] - quadrant->p(0, 7-ii))*
                                          (i.pos[1] - quadrant->p(1, 7-ii))*
                                          (i.pos[2] - quadrant->p(2, 7-ii))) / volume;
                    				
                if (! quadrant->is_hanging (ii))
                  rho_fixed[quadrant->gt (ii)] += i.charge*4.0*pi*weigth / vol_patch[quadrant->gt (ii)];
                else
                  for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                    rho_fixed[quadrant->gparent (jj, ii)] += i.charge*4.0*pi*weigth / 
                                                               (quadrant->num_parents (ii) * vol_patch[quadrant->gt (ii)]);
              }
            //break;
            }
      }
  
  distributed_sparse_matrix A, A_in; 
  A.set_ranges (tmsh.num_owned_nodes ());
  A_in.set_ranges (tmsh.num_owned_nodes ());
  
  distributed_vector  rhs (tmsh.num_owned_nodes (), mpicomm);
  distributed_vector  sigma_free_in (tmsh.num_owned_nodes (), mpicomm); 
  
  bim3a_solution_with_ghosts (tmsh, psi);
  bim3a_advection_diffusion (tmsh, epsilon, psi, A);
  bim3a_advection_diffusion (tmsh, epsilon_in, psi, A_in);
  // A += A_in;
  // bim3a_advection_diffusion (tmsh, epsilon_out, psi, A);
  
  bim3a_solution_with_ghosts (tmsh, ones, replace_op);
  bim3a_reaction (tmsh, reaction, ones, A); 
    
  bim3a_solution_with_ghosts (tmsh, rho_fixed);
  bim3a_solution_with_ghosts (tmsh, sigma_free_in);
  bim3a_rhs (tmsh, const_ones, rho_fixed, rhs);
  bim3a_rhs (tmsh, ones_in, rho_fixed, sigma_free_in);
    
  
  // Set boundary conditions.
  dirichlet_bcs3 bcs;
  if (bc == 1) //hom Dir bc
  {
    for (auto const & ibc : bcells){
    	auto cella = ibc.first;
    	auto lato = ibc.second;
    	bcs.push_back (std::make_tuple (cella, lato, 
    		       [] (double x, double y, double z) {return 0;}));
    }
    bim3a_dirichlet_bc (tmsh, bcs, A, rhs);
  }
  
  A.assemble ();
  A_in.assemble (); 
  rhs.assemble();
  sigma_free_in.assemble ();
  
  tmsh.octbin_export ("rho_0", rho_fixed);
  tmsh.octbin_export ("rhs_0", rhs); 
  
  mumps mumps_solver;
  
  std::vector<double> vals;
  std::vector<int> irow, jcol;

  A.aij (vals, irow, jcol, mumps_solver.get_index_base ());

  mumps_solver.set_lhs_distributed ();
  mumps_solver.set_distributed_lhs_structure (tmsh.num_global_nodes (), irow, jcol);
  mumps_solver.set_distributed_lhs_data (vals);
  mumps_solver.set_rhs_distributed (rhs);

  std::cout << "mumps_solver.analyze () = "
            << mumps_solver.analyze ()
            << std::endl;
  std::cout << "mumps_solver.factorize () = "
            << mumps_solver.factorize ()
            << std::endl;
  std::cout << "mumps_solver.solve () = "
            << mumps_solver.solve ()
            << std::endl;         
  
  distributed_vector phi  = mumps_solver.get_distributed_solution ();
  
  bim3a_solution_with_ghosts (tmsh, phi, replace_op);

  tmsh.octbin_export ("phi_0", phi);
  
  ///////
  
  //phi.assemble ();
  //distributed_vector  tmp_vec (tmsh.num_owned_nodes (), mpicomm); 
  auto tmp_vec=A_in*phi;
  //bim3a_rhs (tmsh, const_ones, rho_fixed, tmp_vec);
  tmp_vec.assemble ();
  
  std::transform(sigma_free_in.get_owned_data ().begin (), 
                 sigma_free_in.get_owned_data ().end (),
                 tmp_vec.get_owned_data ().begin (), 
                 sigma_free_in.get_owned_data ().begin (),
                 [] (double a, double b) {return b-a;});
  
  sigma_free_in.assemble ();
  tmsh.octbin_export ("sigma_0", sigma_free_in);
  
  mumps_solver.cleanup ();
}

*/

void
poisson_boltzmann::lis_compute_electric_potential ()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  if (rank == 0)
     std::cout << "\nStarting LIS solution" << std::endl;
     
  // diffusion
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/(e*e);   //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/(e*e); //adim e_out
  epsilon.assign (tmsh.num_local_quadrants (), eps_in); //e_in
  
  distributed_vector  psi (tmsh.num_owned_nodes ());
  psi.get_owned_data ().assign (psi.get_owned_data ().size (), 0.0);
    
  for (auto epsp = epsilon.begin (), mp = marker.begin ();
       epsp != epsilon.end () || mp != marker.end ();
       ++epsp, ++mp) 
    if ((*mp) < 0.6)
      (*epsp) = eps_in; 
    else  
      (*epsp) = eps_out;
        

  /////////////////////////////////////////////////////////           
  //reactions
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/(e_0*e_out*kb*T);  
  reaction.assign (tmsh.num_local_quadrants (), 0.0);
  
  distributed_vector ones (tmsh.num_owned_nodes ()); 
  ones.get_owned_data ().assign (ones.get_owned_data ().size (), 1.0); 
  
  for (auto rp = reaction.begin (), mp = marker.begin ();
       rp != reaction.end () || mp != marker.end ();
       ++rp, ++mp)
    if ((*mp) > 0.6) 
      (*rp) = eps_out*k2;
      
 
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  rho_fixed = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ()); 
  rho_fixed->get_owned_data ().assign (tmsh.num_owned_nodes (), 0.0); 
  std::vector<double> const_ones (tmsh.num_local_quadrants (), 1.0);
  
  distributed_vector  vol_patch (tmsh.num_owned_nodes (), mpicomm); 
  bim3a_solution_with_ghosts (tmsh, ones, replace_op);
  bim3a_rhs (tmsh, const_ones, ones, vol_patch);
  vol_patch.assemble (); 
  

  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
           quadrant != this->tmsh.end_quadrant_sweep ();
           ++quadrant)
  {
    if (marker[quadrant->get_forest_quad_idx()] < 0.6)
    {  // only inside the molecule                                  
      for (const NS::Atom& i : atoms) 
        if (is_in (i, quadrant)) 
        {
          //linear approx:
          double volume = (quadrant->p(0, 7) - quadrant->p(0, 0)) *
                          (quadrant->p(1, 7) - quadrant->p(1, 0)) *
                          (quadrant->p(2, 7) - quadrant->p(2, 0)); //volume
          for (int ii = 0; ii < 8; ++ii)
          {
            double weigth = std::abs ((i.pos[0] - quadrant->p(0, 7-ii))*
                                      (i.pos[1] - quadrant->p(1, 7-ii))*
                                      (i.pos[2] - quadrant->p(2, 7-ii))) / volume;
                				
            if (! quadrant->is_hanging (ii))
              (*rho_fixed)[quadrant->gt (ii)] += i.charge*4.0*pi*weigth / vol_patch[quadrant->gt (ii)];
            else
              for (int jj = 0; jj < quadrant->num_parents (ii); ++jj)
                (*rho_fixed)[quadrant->gparent (jj, ii)] += i.charge*4.0*pi*weigth / 
                                                           (quadrant->num_parents (ii) * vol_patch[quadrant->gt (ii)]);
            
          }
          //break;
        }
		}   
  }

  
  //////////////////////////////////////////////////////////////////
  MPI_Barrier(mpicomm);
  auto start3 = std::chrono::steady_clock::now();

  distributed_sparse_matrix A; 
  A.set_ranges (tmsh.num_owned_nodes ());
  
  distributed_vector  rhs (tmsh.num_owned_nodes (), mpicomm);
  
  bim3a_solution_with_ghosts (tmsh, psi);
   
  // bim3a_advection_diffusion (tmsh, epsilon, psi, A);
  bim3a_laplacian (tmsh, epsilon, A);

  bim3a_solution_with_ghosts (tmsh, ones, replace_op);
  bim3a_reaction (tmsh, reaction, ones, A); 
    
  bim3a_solution_with_ghosts (tmsh, *rho_fixed);
  bim3a_rhs (tmsh, const_ones, *rho_fixed, rhs);
  
  MPI_Barrier(mpicomm);
  auto end3 = std::chrono::steady_clock::now(); 
  if(rank==0)
  {
    std::cout << "\nTime to calculate elements of A and rhs:  "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end3-start3).count ()
              << " ms"
              <<std::endl;      
  }

  // Set boundary conditions.
  dirichlet_bcs3 bcs;
  if (bc == 1) //hom Dir bc
  {
    MPI_Barrier(mpicomm);
    auto start = std::chrono::steady_clock::now();
    
    for (auto const & ibc : bcells){
    	auto cella = ibc.first;
    	auto lato = ibc.second;
    	bcs.push_back (std::make_tuple (cella, lato, 
    								 [] (double x, double y, double z) {return 0;}));
    }
    bim3a_dirichlet_bc (tmsh, bcs, A, rhs);
    
    MPI_Barrier(mpicomm);
    auto end = std::chrono::steady_clock::now(); 
    if(rank==0)
    {
      std::cout << "\nTime for impose BC:  "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count ()
                << " ms"
                <<std::endl;      
    }
  }
  
  if (bc == 2) //coulombic Dir bc 
  {
    MPI_Barrier(mpicomm);
    auto start = std::chrono::steady_clock::now();
    for (auto const & ibc : bcells){
    	auto cella = ibc.first;
    	auto lato = ibc.second;
    	bcs.push_back (std::make_tuple (cella, lato, 
    		             [&] (double x, double y, double z) {return coulomb_boundary_conditions (x,y,z);}));
    }
    bim3a_dirichlet_bc (tmsh, bcs, A, rhs);
    
    MPI_Barrier(mpicomm);
    auto end = std::chrono::steady_clock::now(); 
    if(rank==0)
    {
      std::cout << "\nTime for impose BC:  "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count ()
                << " ms"
                <<std::endl;      
    }
  }
  /////////////////////////////////////////////////////////

  MPI_Barrier(mpicomm);
  auto start1 = std::chrono::steady_clock::now();
  
  A.assemble ();
  rhs.assemble();
  
  MPI_Barrier(mpicomm);
  auto end1 = std::chrono::steady_clock::now(); 
  if(rank==0)
  {
    std::cout << "\nTime to assemble matrix A and rhs:  "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1-start1).count ()
              << " ms"
              <<std::endl;      
  }
  
  MPI_Barrier(mpicomm);
  auto start_sol = std::chrono::steady_clock::now();
  

  //CSR  
  std::vector<double> vals;
  std::vector<int> irow, jcol;

  
  A.csr(vals, jcol, irow);

  // lis RHS
  LIS_INT i, is, ie, n_rhs, ln; 
  LIS_VECTOR rhs_lis;
  //n_rhs = tmsh.num_global_nodes();
  ln = rhs.get_owned_data().size(); 
  
  lis_vector_create(mpicomm, &rhs_lis);
  lis_vector_set_size(rhs_lis, ln, 0);
  lis_vector_get_range(rhs_lis, &is, &ie);
  
  for (i=is; i<ie; i++)
    lis_vector_set_value (LIS_INS_VALUE, i, rhs.get_owned_data()[i-is], rhs_lis); 
  //lis_vector_print(rhs_lis);
  // lis PHI
  LIS_VECTOR phi_lis;
  
  lis_vector_create(mpicomm, &phi_lis);
  lis_vector_set_size(phi_lis, ln, 0);
  //lis_vector_set_size(phi_lis, 0, n_rhs);
  lis_vector_get_range(phi_lis, &is, &ie);
  
  // lis MATRIX
  LIS_INT n, nnz; //n: matrix dim ; nnz: numb of non zero elems 
  LIS_INT *index; //array of integer containing the col index of non zero elems
  LIS_INT *ptr; //array of integer with starting points of rows 
  LIS_SCALAR *value; //array of double stores non-zero elements of matrix A along the row
  LIS_MATRIX A_lis; //array of integer containing the col index of non zero elems 
  
  nnz = A.owned_nnz();
  n = tmsh.num_owned_nodes ();
  
  lis_matrix_create(mpicomm, &A_lis);
  
  lis_matrix_set_size(A_lis, n, 0);
  ptr = &irow[0];
  index = &jcol[0];
  value = &vals[0];
  
  lis_matrix_set_csr(nnz, ptr, index, value, A_lis);
  
  lis_matrix_assemble(A_lis);
  
  //Solve linear system
  LIS_SOLVER solver;
  
  lis_solver_create(&solver);
  
  std::string opts = linear_solver_options;
  
  lis_solver_set_option(&opts[0], solver);
  
  lis_solve(A_lis, rhs_lis, phi_lis, solver);
  phi = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ());

  lis_vector_get_values(phi_lis, is, ln, phi->get_owned_data ().data());
  
  /*
  distributed_vector phi (tmsh.num_owned_nodes (), mpicomm); 
  {
    lis_distributed solver;
    A.csr (vals, jcol, irow, solver.get_index_base ());
    solver.set_lhs_structure (tmsh.num_global_nodes (), irow, jcol);
    solver.set_lhs_data (vals);
    solver.set_rhs (rhs.get_owned_data ());
    solver.set_initial_guess (phi.get_owned_data ());
    solver.set_iterative_method ("Conjugate Gradient");
    solver.set_preconditioner ("ilu");
    solver.set_convergence_condition ("nrm1_b");
    solver.analyze ();
    solver.factorize ();
    solver.solve ();
    solver.cleanup ();
  }
  */
  MPI_Barrier(mpicomm);
  auto end_sol = std::chrono::steady_clock::now(); 
  if(rank==0)
  {
    std::cout << "\nTime to solve linear problem:  "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end_sol- start_sol).count ()
              << " ms"
              <<std::endl;      
  }

  bim3a_solution_with_ghosts (tmsh, *phi, replace_op);
  phi_energy = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ());
  phi_energy -> get_owned_data () = phi->get_owned_data();

  rho_fixed_energy = std::make_unique<distributed_vector> (tmsh.num_owned_nodes ());
  rho_fixed_energy -> get_owned_data () = rho_fixed->get_owned_data();
  ///////
  
	/////////////////////////////////////////////////////////
	
//   MPI_Barrier(mpicomm);
//   auto start4 = std::chrono::steady_clock::now(); 
//   tmsh.octbin_export_quadrant ("epsilon_0", epsilon);
//   tmsh.octbin_export ("phi_0", phi_tmp);
//   tmsh.octbin_export ("rho_0", rho_fixed);
//   MPI_Barrier(mpicomm);
//   auto end4 = std::chrono::steady_clock::now(); 
//   if(rank==0)
//   {
//     std::cout << "\nTime to export files:  "
//               << std::chrono::duration_cast<std::chrono::milliseconds>(end4-start4).count ()
//               << " ms"
//               <<std::endl;      
//   }
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////



void 
poisson_boltzmann::surface_integrals_energy()
{
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  if (rank == 0)
     std::cout << "\nStarting energy calculation with surface integrals" << std::endl;
  
  
  double net_charge = 0.0;
  for (const NS::Atom& i : atoms){
    net_charge += i.charge;
  }

  ////////////////////////////////////////////////////////
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/(e*e);   //adim e_in
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/(e*e); //adim e_out
  epsilon_in.assign (tmsh.num_local_quadrants (), 0.0);
  ones_in.assign (tmsh.num_local_quadrants (), 0.0); 

    
  for (auto epsp_in = epsilon_in.begin (),  
            onesp_in = ones_in.begin (), 
            mp = marker.begin ();
       epsp_in != epsilon_in.end () 
       || onesp_in != ones_in.end ()
       || mp != marker.end ();
       ++epsp_in, ++onesp_in, ++mp) 
    if ((*mp) < 0.6)
    {
      (*epsp_in) = eps_in;
      (*onesp_in) = 1.0; 
    }
    

  distributed_vector tmp_sf    (tmsh.num_owned_nodes (), mpicomm); 
  distributed_vector sf_nodes  (tmsh.num_owned_nodes (), mpicomm); // surface nodes
  distributed_vector coord_X   (tmsh.num_owned_nodes (), mpicomm);
  distributed_vector coord_Y   (tmsh.num_owned_nodes (), mpicomm);
  distributed_vector coord_Z   (tmsh.num_owned_nodes (), mpicomm);

  tmp_sf.get_owned_data ().assign (tmp_sf.get_owned_data ().size (), 0.0);
  sf_nodes.get_owned_data ().assign (sf_nodes.get_owned_data ().size (), 0.0);

  coord_X.get_owned_data ().assign (coord_X.get_owned_data ().size (), 0.0);
  coord_Y.get_owned_data ().assign (coord_Y.get_owned_data ().size (), 0.0);
  coord_Z.get_owned_data ().assign (coord_Z.get_owned_data ().size (), 0.0);
  

  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
           quadrant != this->tmsh.end_quadrant_sweep ();
           ++quadrant)
  {
    for (int kk = 0; kk < 8; ++kk)
    {
      if (! quadrant->is_hanging (kk))
      {
        coord_X[quadrant->gt (kk)] = quadrant->p(0,kk);
        coord_Y[quadrant->gt (kk)] = quadrant->p(1,kk);
        coord_Z[quadrant->gt (kk)] = quadrant->p(2,kk);
      }
    }
    

    if (marker[quadrant->get_forest_quad_idx()] < 0.6)
    {  // only inside the molecule                         
      for (int kk = 0; kk < 8; ++kk)
      {
        if (! quadrant->is_hanging (kk))
        {
          tmp_sf[quadrant->gt (kk)] += 1;
        }
      }
    }
  }
  tmp_sf.assemble();
  coord_X.assemble(replace_op);
  coord_Y.assemble(replace_op);
  coord_Z.assemble(replace_op);

  for (auto jj = 0; jj< tmsh.num_owned_nodes (); ++jj){
    if (tmp_sf.get_owned_data ()[jj]==0 || tmp_sf.get_owned_data ()[jj] == 8)
    {
      sf_nodes.get_owned_data ()[jj] = 0.0;
    } else {
      sf_nodes.get_owned_data ()[jj] = 1.0;
    }
  }
  // sf_nodes.assemble(replace_op);
  
  bim3a_solution_with_ghosts (tmsh, sf_nodes, replace_op);
  tmsh.octbin_export ("sf_nodes_0", sf_nodes);

  distributed_vector  psi (tmsh.num_owned_nodes ());
  psi.get_owned_data ().assign (psi.get_owned_data ().size (), 0.0);
  bim3a_solution_with_ghosts (tmsh, psi);
  bim3a_solution_with_ghosts (tmsh, *phi_energy,replace_op);
  bim3a_solution_with_ghosts (tmsh, *rho_fixed_energy, replace_op);
  distributed_sparse_matrix A_in; 
  A_in.set_ranges (tmsh.num_owned_nodes ());
  bim3a_advection_diffusion (tmsh, epsilon_in, psi, A_in);
  A_in.assemble ();
  
  distributed_vector  sigma_free_in (tmsh.num_owned_nodes (), mpicomm); 
  sigma_free_in.get_owned_data().assign(tmsh.num_owned_nodes (), 0.0);
  bim3a_rhs (tmsh, ones_in, *rho_fixed_energy, sigma_free_in);
  sigma_free_in.assemble();

  distributed_vector tmp_vec (tmsh.num_owned_nodes (), mpicomm);

  tmp_vec = A_in*(*phi_energy);
  tmp_vec.assemble ();
  ////////////////////////////////////////////////////////////////////////////
  
  // for (auto iA = (*phi_energy).get_range_start (); iA < (*phi_energy).get_range_end (); ++iA)
  //     for (auto j = A_in[iA].begin (); j != A_in[iA].end (); ++j)
  //       if (A_in.col_val (j) != 0)
  //         tmp_vec(iA) += A_in.col_val (j) * (*phi_energy)(A_in.col_idx (j));
  
  // tmp_vec.assemble (replace_op);
  ////////////////////////////////////////////////////////////////////////////
             
  std::transform(sigma_free_in.get_owned_data ().begin (), 
                 sigma_free_in.get_owned_data ().end (),
                 tmp_vec.get_owned_data ().begin (), 
                 sigma_free_in.get_owned_data ().begin (),
                 [] (double a, double b) {return a-b;});

  // sigma_free_in.assemble(replace_op);

  bim3a_solution_with_ghosts (tmsh, sigma_free_in, replace_op);
  tmsh.octbin_export ("sigma_0", sigma_free_in);
  tmsh.octbin_export ("phi_en_0", *phi_energy);
  
  double charge_pol = std::accumulate(sigma_free_in.get_owned_data ().begin (), 
                                      sigma_free_in.get_owned_data ().end (),
                                      0.0);

  double energy_pol    = 0.0;
  double energy_react1 = 0.0;
  double energy_react  = 0.0;
  double distance      = 0.0;
  double first_int     = 0.0;
  double second_int    = 0.0;

  double coul_energy   = 0.0;
  double den_in        = 1.0/(eps_in);
  int i_atom = 0;
  int j_atom = 0;

  
  distributed_vector u_i (tmsh.num_owned_nodes (), mpicomm);
  distributed_vector v_i (tmsh.num_owned_nodes (), mpicomm);
  bim3a_solution_with_ghosts (tmsh, u_i, replace_op); 

  const double constant = 0.5*(1.0/eps_out - 1.0/eps_in)/(4.0*pi);
  const double constant_react = 0.5/(4.0*pi*eps_out);


  for (const NS::Atom& i : atoms){
    if ( std::fabs(i.charge) > 0.0){
      u_i.get_owned_data ().assign (u_i.get_owned_data ().size (), 0.0); 
      
      for (auto jj = 0; jj< tmsh.num_owned_nodes (); ++jj){
        distance = std::hypot ((i.pos[0] - coord_X.get_owned_data ()[jj]), 
                               (i.pos[1] - coord_Y.get_owned_data ()[jj]), 
                               (i.pos[2] - coord_Z.get_owned_data ()[jj])); 
        u_i.get_owned_data ()[jj] = - 1.0/(eps_in*distance);
        if(sf_nodes.get_owned_data ()[jj] == 1.0)
        { 
          first_int +=  i.charge*sigma_free_in.get_owned_data ()[jj]/distance; 
        }      
      }
      u_i.assemble (replace_op);  

      ///////////////////////////////////////////
      
      v_i=A_in*u_i;  
      v_i.assemble ();

      // for (auto iA = (*phi_energy).get_range_start (); iA < (*phi_energy).get_range_end (); ++iA)
      //   for (auto j = A_in[iA].begin (); j != A_in[iA].end (); ++j)
      //     if (A_in.col_val (j) != 0)
      //       v_i(iA) += A_in.col_val (j) * u_i(A_in.col_idx (j));
   
      // v_i.assemble (replace_op);
      
      // ///////////////////////////////////////////
      
      for (auto jj = 0; jj< tmsh.num_owned_nodes (); ++jj){
        if(sf_nodes.get_owned_data ()[jj] == 1.0){
          second_int += 0.5*i.charge*(*phi_energy).get_owned_data ()[jj]*
                          v_i.get_owned_data ()[jj]/(4.0*pi);                
        }       
      }

    }
  }
  std::cout<<second_int << std::endl;
  energy_pol = constant*first_int;
  energy_react = second_int - constant_react*first_int;
  coul_energy = coul_energy*den_in;
  ///////////////////////////////////////////////////////// 
  
  /////////////////////////////////////////////////////////
  //
  //coulombic energy

  if (rank == 0) {
    for (const NS::Atom& i : atoms){
      if ( std::fabs(i.charge) > 0.0){
        for (const NS::Atom& j : atoms){
          if ( std::fabs(j.charge) > 0.0){ 
            if(j_atom > i_atom ){
              distance = std::hypot((i.pos[0] - j.pos[0]), 
                                  (i.pos[1] - j.pos[1]), 
                                  (i.pos[2] - j.pos[2]));
              coul_energy += i.charge*j.charge/distance;                           
            }
          }
          j_atom ++;
        }
      }
      i_atom ++;
      j_atom = 0;
    }
    
    coul_energy = coul_energy*den_in;
  }  
    

  if (rank == 0) {
  MPI_Reduce(MPI_IN_PLACE, &energy_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
             mpicomm);
  MPI_Reduce(MPI_IN_PLACE, &energy_react, 1, MPI_DOUBLE, MPI_SUM, 0,
             mpicomm);
  MPI_Reduce(MPI_IN_PLACE, &charge_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
             mpicomm);
  } else {
  MPI_Reduce(&energy_pol, &energy_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
             mpicomm);
  MPI_Reduce(&energy_react, &energy_react, 1, MPI_DOUBLE, MPI_SUM, 0,
             mpicomm);
  MPI_Reduce(&charge_pol, &charge_pol, 1, MPI_DOUBLE, MPI_SUM, 0,
             mpicomm);
  }
              
  // Print the result
  if (rank == 0) {
    std::cout << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;   
    std::cout << "Net charge: "
              << std::setprecision(16)<<net_charge
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;   
    std::cout << "Polarization charge: "
              << std::setprecision(16)<<charge_pol/(4.0*pi)
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;   
    std::cout << "Polarization energy value: "
              << std::setprecision(16)<<energy_pol
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;
  
    std::cout << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;   
    std::cout << "Reaction energy value: "
              << std::setprecision(16)<<energy_react 
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;   
    std::cout << "Coulumbic energy value: "
              << std::setprecision(16)<<coul_energy
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;   
    std::cout << "Total energy value: "
              << std::setprecision(16)<<energy_react + coul_energy + energy_pol
              << std::endl;
    std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;
  }
}

double
poisson_boltzmann::coulomb_boundary_conditions(double x, double y, double z)
{
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/(e*e); //adim e_out
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/(e_0*e_out*kb*T); 
  
  double dist = 0.0;
  double pot  = 0.0;
  double k    = std::sqrt(k2);
  for (const NS::Atom& i : atoms)
  {
    dist = std::hypot ((i.pos[0] - x), (i.pos[1] - y), (i.pos[2] - z));
    pot += i.charge*exp(-k*dist)/(dist*eps_out);
  }
  return pot;
}


void 
poisson_boltzmann::grid_energy_electric_field(distributed_vector &phi, std::vector<double> &epsilon, std::string int_rule, double &energy)
{
	double energy_partition = 0.0;
  
  int rank;
  MPI_Comm_rank (mpicomm, &rank);
  if (rank == 0)
     std::cout << "\nStarting grid energy solution" << std::endl;
  if(int_rule == "gauss")
		{
			for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
		         quadrant != this->tmsh.end_quadrant_sweep ();
		         ++quadrant)
		  //sum on all quadrant
				{
					double
						x[2] = {quadrant->p (0,0), quadrant->p (0,7)},
						y[2] = {quadrant->p (1,0), quadrant->p (1,7)},
						z[2] = {quadrant->p (2,0), quadrant->p (2,7)};
					
					double phi_loc[8] = {0,0,0,0,0,0,0,0};
					double eps_loc = epsilon[quadrant->get_forest_quad_idx()];
					for (int ii = 0; ii < 8; ++ii)
						{
							if (! quadrant->is_hanging (ii)){
				  			phi_loc[ii] = phi[quadrant ->gt(ii)];
				  		}
							else{
						    int np = quadrant->num_parents(ii);
						  	for (int pp = 0; pp < np; ++pp)
						  		{
						  			phi_loc[ii] += phi[quadrant ->gparent (pp, ii)];
						   		}
						  	phi_loc[ii] /= np;
							}
						}
				  
					auto fun =
						[&x, &y, &z, &eps_loc, &phi_loc]
						(double X, double Y, double Z) -> double
						  { return 0.5*eps_loc*
						           (dudx(X, Y, Z, x, y, z, phi_loc)*dudx(X, Y, Z, x, y, z, phi_loc) +
						            dudy(X, Y, Z, x, y, z, phi_loc)*dudy(X, Y, Z, x, y, z, phi_loc) +
						            dudz(X, Y, Z, x, y, z, phi_loc)*dudz(X, Y, Z, x, y, z, phi_loc));};
				  energy_partition += quad_integral (x, y, z, fun)/(4.0*pi);
				  
				}
				
		}

	else std::cout << "Errore" << std::endl;
	
	
	MPI_Reduce(&energy_partition, &energy, 1, MPI_DOUBLE, MPI_SUM, 0,
             mpicomm);
             
	
	// Print the result
  if (rank == 0) {
    std::cout << std::endl;
  	std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;   
  	std::cout << "Total energy value: "
   	 	        << std::setprecision(16)<<energy
    	        << std::endl;
  	std::cout <<"+++++++++++++++++++++++++++++++++" << std::endl;
  	std::cout << std::endl;
  }
  ////////////////////////////////////////////////////////////////
}



void 
poisson_boltzmann::analitic_potential()
{
  double eps_out = 4.0*pi*e_0*e_out*kb*T*Angs/(e*e); //adim e_out
  double eps_in = 4.0*pi*e_0*e_in*kb*T*Angs/(e*e);  //adim e_in
  double C_0 = 1.0e3*N_av*ionic_strength; //Bulk concentration of monovalent species
  double k2 = 2.0*C_0*Angs*Angs*e*e/(e_0*e_out*kb*T);    
  double k    = std::sqrt(k2);
  double pippo_out = eps_out*(1.0+k*10.0);
  double pippo_in  = eps_in*10.0; 
  double dist= 0.0;
  
  distributed_vector phi_an (tmsh.num_owned_nodes (), mpicomm);
  distributed_vector field_an (tmsh.num_owned_nodes (), mpicomm);
  bim3a_solution_with_ghosts (tmsh, phi_an, replace_op);
  bim3a_solution_with_ghosts (tmsh, field_an, replace_op);
  
  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
           quadrant != this->tmsh.end_quadrant_sweep ();
           ++quadrant){            
    for (const NS::Atom& i : atoms){
      if ( std::fabs(i.charge) > std::numeric_limits<double>::epsilon () ){
				for (int ii = 0; ii < 8; ++ii){
				  if (! quadrant->is_hanging (ii)){
				    dist = std::hypot ((i.pos[0] - quadrant->p(0,ii)), 
						                   (i.pos[1] - quadrant->p(1,ii)), 
						                   (i.pos[2] - quadrant->p(2,ii))); 
						if(dist <= 10.0){
						  phi_an[quadrant->gt (ii)] = 1.0/(pippo_out*10.0) + 
				                                  i.charge*(10.0-dist)/(pippo_in*dist);
				      field_an[quadrant->gt (ii)] = (1.0 + (10.0-dist)/dist)*
				                                     i.charge/(pippo_in*dist);
						}
						else
				      phi_an[quadrant->gt (ii)] = i.charge*exp(k*(10.0-dist))/(dist*pippo_out);		
				      field_an[quadrant->gt (ii)] = phi_an[quadrant->gt (ii)]*(k +1.0/dist);		    
				  }
				  else{
				    for (int jj = 0; jj < quadrant->num_parents (ii); ++jj){
				      phi_an[quadrant->gparent (jj, ii)] += 0.0;
				      field_an[quadrant->gparent (jj, ii)] += 0.0;
				    } 
				  }
				} 
	    }
	  }        
  }

  phi_an.assemble (replace_op);
  tmsh.octbin_export ("phi_an_0", phi_an);
  field_an.assemble (replace_op);
  tmsh.octbin_export ("field_an_0", field_an);
}


void 
poisson_boltzmann::abs_value_field(distributed_vector &phi)
{
  distributed_vector abs_field (tmsh.num_owned_nodes (), mpicomm);
  bim3a_solution_with_ghosts (tmsh, abs_field, replace_op); 
  
  
  for (auto quadrant = this->tmsh.begin_quadrant_sweep ();
           quadrant != this->tmsh.end_quadrant_sweep ();
           ++quadrant){
    double
      x[2] = {quadrant->p (0,0), quadrant->p (0,7)},
      y[2] = {quadrant->p (1,0), quadrant->p (1,7)},
      z[2] = {quadrant->p (2,0), quadrant->p (2,7)}; 
      double phi_loc[8] = {0,0,0,0,0,0,0,0};
    for(int ii=0; ii<8 ; ++ii){
      if (! quadrant->is_hanging (ii)){
			  phi_loc[ii] = phi[quadrant ->gt(ii)];
      }
			else{
			  int np = quadrant->num_parents(ii);
				for (int pp = 0; pp < np; ++pp)
				{
					phi_loc[ii] += phi[quadrant ->gparent (pp, ii)];
	   		}
		  	phi_loc[ii] /= np;
		  }
      abs_field[quadrant ->gt(ii)] = 
        std::hypot(dudx(x[0], y[0], z[0], x, y, z, phi_loc),
                   dudy(x[0], y[0], z[0], x, y, z, phi_loc),
                   dudz(x[0], y[0], z[0], x, y, z, phi_loc));
    }
    
    
  
  }
  abs_field.assemble (replace_op);
  tmsh.octbin_export ("field_0", abs_field);
  
  
}




