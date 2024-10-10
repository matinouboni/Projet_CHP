#include <iostream>
#include "TimeScheme.h"
#include <fstream>
#include <string>

int main(int argc, char** argv)
{

   if (argc < 2)
   {
      printf("Veuillez entrer le fichier data (exemple : ../Data/data_cas_1.toml).\n");
      exit(0);
   }

   MPI_Init(&argc, &argv);

      // t1 et t2 pour le temps de calcul
      double t1, t2;

      const std::string data_file_name = argv[1];

      DataFile* data_file = new DataFile(data_file_name);
      Function* fct = new Function(data_file);
      Advection* adv = new Advection(fct, data_file);
      TimeScheme* time_scheme = NULL;

      int int_euler_scheme(data_file->Get_time_scheme());
      // Switch pour le choix du schema
      switch(int_euler_scheme)
      {
         case 1:
         time_scheme = new ExplicitScheme(data_file, adv);
         break;
         case 2:
         time_scheme = new ImplicitScheme(data_file, adv);
         break;
         default:
         printf("Ce choix n'est pas possible ! Veuillez recommencer !");
         exit(0);
      }

      int Np, Me, iBeg, iEnd;
      MPI_Status status;
      MPI_Comm_rank(MPI_COMM_WORLD, &Me);
      MPI_Comm_size(MPI_COMM_WORLD, &Np);
      
      // Repartion de la charge
      charge(Me, data_file->Get_Ny()*data_file->Get_Nx(), Np, &iBeg, &iEnd);

      // Debut du chrono
      t1 = MPI_Wtime();

      // Vecteur solution local (pour un proc Me) 
      std::vector<double> U(iEnd-iBeg+1,0.0);

      // Initialisation de la solution
      adv->InitialCondition(U, iBeg, iEnd);
      double t(0.0);
      int it(0);

      // Sauvegarde de la solution numerique et exacte
      time_scheme->SaveSol(U, "sol", it, iBeg, iEnd, Me);
      time_scheme->SaveSol(adv->ExactSol(t, iBeg, iEnd), "ex", it, iBeg, iEnd, Me);
      it++;

      // Boucle en temps
      while (t < data_file->Get_Tf()){

         // Calcul de la solution en tn+1
         time_scheme->Integrate(t, U, iBeg, iEnd, Me, Np, status);

         // Sauvgarde des solutions
         time_scheme->SaveSol(U, "sol", it, iBeg, iEnd, Me);
         time_scheme->SaveSol(adv->ExactSol(t, iBeg, iEnd), "ex", it, iBeg, iEnd, Me);
         it++;
      }

      // Fin du chrono
      t2 = MPI_Wtime();

      printf("Temps proc %d = %lf secondes\n", Me, t2-t1);

      delete data_file, delete fct, delete adv, delete time_scheme;

   MPI_Finalize();

   return 0;
}
