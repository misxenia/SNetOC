% publish_all demo files in html

opts.format = 'html';
opts.outputDir = 'docs';

publish('demo_overlappingcommunity', opts);
publish('demo_polblogs', opts);
publish('demo_sparsity', opts);
publish('demo_simulations', opts);
publish('demo_usairport', opts);
